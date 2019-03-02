/******************************************************************************/
/*                       LYNX, a MOOSE-based application                      */
/*                                                                            */
/*          Copyright (C) 2017 by Antoine B. Jacquey and Mauro Cacace         */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*    This program is free software: you can redistribute it and/or modify    */
/*    it under the terms of the GNU General Public License as published by    */
/*      the Free Software Foundation, either version 3 of the License, or     */
/*                     (at your option) any later version.                    */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program. If not, see <http://www.gnu.org/licenses/>     */
/******************************************************************************/

#include "LynxSolidMomentum.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "Assembly.h"
#include "SystemBase.h"
#include "libmesh/quadrature.h"

registerMooseObject("LynxApp", LynxSolidMomentum);

template <>
InputParameters
validParams<LynxSolidMomentum>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Solid momentum kernel.");
  params.addRequiredCoupledVar(
      "displacements", "The string of displacements variables suitable for the problem statement.");
  params.addCoupledVar("temperature", "The temperature variable.");
  params.addCoupledVar("fluid_pressure", "The fluid pressure variable.");
  params.addCoupledVar("lithostatic_pressure", "The lithostatic pressure variable.");
  params.addCoupledVar("dynamic_pressure", "The dynamic pressure variable.");
  params.set<bool>("use_displaced_mesh") = false;
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in (0 for x, "
                                        "1 for y, 2 for z).");
  params.addParam<bool>(
      "volumetric_locking_correction", false, "Flag to correct volumetric locking");
  return params;
}

LynxSolidMomentum::LynxSolidMomentum(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _ndisp(coupledComponents("displacements")),
    _coupled_temp(isCoupled("temperature")),
    _coupled_pf(isCoupled("fluid_pressure")),
    _pf(_coupled_pf ? coupledValue("fluid_pressure") : _zero),
    _coupled_plith(isCoupled("lithostatic_pressure")),
    _coupled_pdyn(isCoupled("dynamic_pressure")),
    _component(getParam<unsigned int>("component")),
    _vol_locking_correction(getParam<bool>("volumetric_locking_correction")),
    _stress(getMaterialProperty<RankTwoTensor>("stress")),
    _biot(getDefaultMaterialProperty<Real>("biot_coefficient")),
    _tangent_modulus(getMaterialProperty<RankFourTensor>("tangent_modulus")),
    _dthermal_strain_dtemp(getDefaultMaterialProperty<RankTwoTensor>("dthermal_strain_dtemp")),
    _gravity(getDefaultMaterialProperty<RealVectorValue>("gravity_vector")),
    _rho_b(getDefaultMaterialProperty<Real>("bulk_density")),
    _drho_dtemp(getDefaultMaterialProperty<Real>("drho_dtemp")),
    _drho_dev(getDefaultMaterialProperty<Real>("drho_dev")),
    _avg_grad_test(_test.size(), std::vector<Real>(3, 0.0)),
    _avg_grad_phi(_phi.size(), std::vector<Real>(3, 0.0)),
    _disp_var(_ndisp),
    _temp_var(_coupled_temp ? coupled("temperature") : 0),
    _pf_var(_coupled_pf ? coupled("fluid_pressure") : 0),
    _plith_var(_coupled_plith ? coupled("lithostatic_pressure") : 0),
    _pdyn_var(_coupled_pdyn ? coupled("dynamic_pressure") : 0)
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError("The number of displacement variables supplied must match the mesh dimension.");

  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

void
LynxSolidMomentum::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  _local_re.resize(re.size());
  _local_re.zero();

  if (_vol_locking_correction)
    computeAverageGradientTest();

  precalculateResidual();
  for (_i = 0; _i < _test.size(); ++_i)
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      _local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual();

  re += _local_re;

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (const auto & var : _save_in)
      var->sys().solution().add_vector(_local_re, var->dofIndices());
  }
}

Real
LynxSolidMomentum::computeQpResidual()
{
  RealVectorValue stress_row = _stress[_qp].row(_component);
  stress_row(_component) += _biot[_qp] * _pf[_qp];
  RealVectorValue grav_term = -_rho_b[_qp] * _gravity[_qp];

  Real residual = stress_row * _grad_test[_i][_qp] + grav_term(_component) * _test[_i][_qp];

  if (_vol_locking_correction)
    residual += (_stress[_qp].trace() / 3.0 + _biot[_qp] * _pf[_qp]) *
                (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component));

  return residual;
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

void
LynxSolidMomentum::computeJacobian()
{
  if (_vol_locking_correction)
  {
    computeAverageGradientTest();
    computeAverageGradientPhi();
  }

  Kernel::computeJacobian();
}

Real
LynxSolidMomentum::computeQpJacobian()
{
  Real jacobian = 0.0;
  jacobian += elasticJacobian(
      _tangent_modulus[_qp], _component, _component, _grad_test[_i][_qp], _grad_phi[_j][_qp]);

  RealVectorValue dgrav_term = -_drho_dev[_qp] * _gravity[_qp];
  jacobian += dgrav_term(_component) * _test[_i][_qp] * _grad_phi[_j][_qp](_component);

  if (_vol_locking_correction)
  {
    Real sum_C3x3 = _tangent_modulus[_qp].sum3x3();
    RealGradient sum_C3x1 = _tangent_modulus[_qp].sum3x1();
    // jacobian = Bbar^T_i * C * Bbar_j where Bbar = B + Bvol
    // jacobian = B^T_i * C * B_j + Bvol^T_i * C * Bvol_j +  Bvol^T_i * C * B_j + B^T_i * C * Bvol_j

    // Bvol^T_i * C * Bvol_j
    jacobian += sum_C3x3 * (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) *
                (_avg_grad_phi[_j][_component] - _grad_phi[_j][_qp](_component)) / 9.0;

    // B^T_i * C * Bvol_j
    jacobian += sum_C3x1(_component) * _grad_test[_i][_qp](_component) *
                (_avg_grad_phi[_j][_component] - _grad_phi[_j][_qp](_component)) / 3.0;

    // Bvol^T_i * C * B_j
    RankTwoTensor phi;
    if (_component == 0)
    {
      phi(0, 0) = _grad_phi[_j][_qp](0);
      phi(0, 1) = phi(1, 0) = _grad_phi[_j][_qp](1);
      phi(0, 2) = phi(2, 0) = _grad_phi[_j][_qp](2);
    }
    else if (_component == 1)
    {
      phi(1, 1) = _grad_phi[_j][_qp](1);
      phi(0, 1) = phi(1, 0) = _grad_phi[_j][_qp](0);
      phi(1, 2) = phi(2, 1) = _grad_phi[_j][_qp](2);
    }
    else if (_component == 2)
    {
      phi(2, 2) = _grad_phi[_j][_qp](2);
      phi(0, 2) = phi(2, 0) = _grad_phi[_j][_qp](0);
      phi(1, 2) = phi(2, 1) = _grad_phi[_j][_qp](1);
    }

    jacobian += (_tangent_modulus[_qp] * phi).trace() *
                (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) / 3.0;
  }

  return jacobian;
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

void
LynxSolidMomentum::computeOffDiagJacobian(MooseVariableFEBase & jvar)
{
  if (_vol_locking_correction)
  {
    computeAverageGradientPhi();
    computeAverageGradientTest();
  }

  Kernel::computeOffDiagJacobian(jvar);
}

Real
LynxSolidMomentum::computeQpOffDiagJacobian(unsigned int jvar)
{
  // // off-diagonal Jacobian with respect to a coupled displacement component
  for (unsigned int coupled_component = 0; coupled_component < _ndisp; ++coupled_component)
    if (jvar == _disp_var[coupled_component])
    {
      Real jacobian = 0.0;
      jacobian += elasticJacobian(_tangent_modulus[_qp],
                                  _component,
                                  coupled_component,
                                  _grad_test[_i][_qp],
                                  _grad_phi[_j][_qp]);

      RealVectorValue dgrav_term = -_drho_dev[_qp] * _gravity[_qp];
      jacobian += dgrav_term(_component) * _test[_i][_qp] * _grad_phi[_j][_qp](coupled_component);

      if (_vol_locking_correction)
      {
        const Real sum_C3x3 = _tangent_modulus[_qp].sum3x3();
        const RealGradient sum_C3x1 = _tangent_modulus[_qp].sum3x1();

        // jacobian = Bbar^T_i * C * Bbar_j where Bbar = B + Bvol
        // jacobian = B^T_i * C * B_j + Bvol^T_i * C * Bvol_j +  Bvol^T_i * C * B_j + B^T_i * C *
        // Bvol_j

        // Bvol^T_i * C * Bvol_j
        jacobian += sum_C3x3 * (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) *
                    (_avg_grad_phi[_j][coupled_component] - _grad_phi[_j][_qp](coupled_component)) /
                    9.0;

        // B^T_i * C * Bvol_j
        jacobian += sum_C3x1(_component) * _grad_test[_i][_qp](_component) *
                    (_avg_grad_phi[_j][coupled_component] - _grad_phi[_j][_qp](coupled_component)) /
                    3.0;

        // Bvol^T_i * C * B_i
        RankTwoTensor phi;
        for (unsigned int i = 0; i < 3; ++i)
          phi(coupled_component, i) = _grad_phi[_j][_qp](i);

        jacobian += (_tangent_modulus[_qp] * phi).trace() *
                    (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component)) / 3.0;
      }

      return jacobian;
    }

  // if (_coupled_temp && jvar == _temp_var)
  //   return -(_elasticity_tensor[_qp] * _dthermal_strain_dtemp[_qp] *
  //            _grad_test[_i][_qp])(_component)*_phi[_j][_qp];

  if (_coupled_pf && jvar == _pf_var)
    return -_biot[_qp] * _phi[_j][_qp] * _grad_test[_i][_qp](_component);

  if (_coupled_plith && jvar == _plith_var)
    return -_phi[_j][_qp] * _grad_test[_i][_qp](_component);

  if (_coupled_pdyn && jvar == _pdyn_var)
    return -_phi[_j][_qp] * _grad_test[_i][_qp](_component);

  return 0.0;
}

Real
LynxSolidMomentum::elasticJacobian(const RankFourTensor & jacobian_r4t,
                                   unsigned int i,
                                   unsigned int k,
                                   const RealGradient & grad_test,
                                   const RealGradient & grad_phi)
{
  Real sum = 0.0;
  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    for (unsigned int l = 0; l < LIBMESH_DIM; ++l)
      sum += jacobian_r4t(i, j, k, l) * grad_phi(l) * grad_test(j);
  return sum;
}

void
LynxSolidMomentum::computeAverageGradientTest()
{
  // Calculate volume averaged value of shape function derivative
  _avg_grad_test.resize(_test.size());
  for (_i = 0; _i < _test.size(); ++_i)
  {
    _avg_grad_test[_i].resize(3);
    _avg_grad_test[_i][_component] = 0.0;
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      _avg_grad_test[_i][_component] += _grad_test[_i][_qp](_component) * _JxW[_qp] * _coord[_qp];

    _avg_grad_test[_i][_component] /= _current_elem_volume;
  }
}

void
LynxSolidMomentum::computeAverageGradientPhi()
{
  // Calculate volume average derivatives for phi
  _avg_grad_phi.resize(_phi.size());
  for (_i = 0; _i < _phi.size(); ++_i)
  {
    _avg_grad_phi[_i].resize(3);
    for (unsigned int component = 0; component < _mesh.dimension(); ++component)
    {
      _avg_grad_phi[_i][component] = 0.0;
      for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
        _avg_grad_phi[_i][component] += _grad_phi[_i][_qp](component) * _JxW[_qp] * _coord[_qp];

      _avg_grad_phi[_i][component] /= _current_elem_volume;
    }
  }
}