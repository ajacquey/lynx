/******************************************************************************/
/*                            This file is part of                            */
/*                       LYNX, a MOOSE-based application                      */
/*                    Lithosphere dYnamic Numerical toolboX                   */
/*                                                                            */
/*          Copyright (C) 2017 by Antoine B. Jacquey and Mauro Cacace         */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*                Licensed under GNU General Public License 3,                */
/*                       please see LICENSE for details                       */
/*                  or http://www.gnu.org/licenses/gpl.html                   */
/******************************************************************************/

#include "LynxADSolidMomentum.h"
#include "libmesh/quadrature.h"

registerMooseObject("LynxApp", LynxADSolidMomentum);

InputParameters
LynxADSolidMomentum::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Solid momentum kernel.");
  params.addCoupledVar("fluid_pressure", 0.0, "The fluid pressure variable.");
  params.set<bool>("use_displaced_mesh") = false;
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in (0 for x, "
                                        "1 for y, 2 for z).");
  params.addParam<bool>(
      "volumetric_locking_correction", false, "Flag to correct volumetric locking");
  return params;
}

LynxADSolidMomentum::LynxADSolidMomentum(const InputParameters & parameters)
  : ADKernel(parameters),
    _pf(adCoupledValue("fluid_pressure")),
    _component(getParam<unsigned int>("component")),
    _vol_locking_correction(getParam<bool>("volumetric_locking_correction")),
    _stress(getADMaterialProperty<RankTwoTensor>("stress")),
    _coupled_pf(hasADMaterialProperty<Real>("biot_coefficient")),
    _biot(_coupled_pf ? &getADMaterialProperty<Real>("biot_coefficient") : nullptr),
    _coupled_grav(hasMaterialProperty<RealVectorValue>("gravity_vector")),
    _gravity(_coupled_grav ? &getMaterialProperty<RealVectorValue>("gravity_vector") : nullptr),
    _rho_b(_coupled_grav ? &getADMaterialProperty<Real>("bulk_density") : nullptr),
    _avg_grad_test()
{
}

ADReal
LynxADSolidMomentum::computeQpResidual()
{
  ADRealVectorValue stress_row = _stress[_qp].row(_component);
  if (_coupled_pf)
    stress_row(_component) -= (*_biot)[_qp] * _pf[_qp];
  ADRealVectorValue grav_term = _coupled_grav ? -(*_rho_b)[_qp] * (*_gravity)[_qp] : ADRealVectorValue();

  ADReal residual = stress_row * _grad_test[_i][_qp] + grav_term(_component) * _test[_i][_qp];

  // volumetric locking correction
  if (_vol_locking_correction)
    residual += (_avg_grad_test[_i] - _grad_test[_i][_qp](_component)) / 3.0 * _stress[_qp].trace();

  return residual;
}

void
LynxADSolidMomentum::precalculateResidual()
{
  if (!_vol_locking_correction)
    return;

  ADReal ad_current_elem_volume = 0.0;
  for (unsigned int qp = 0; qp < _qrule->n_points(); qp++)
    ad_current_elem_volume += _ad_JxW[qp] * _ad_coord[qp];

  // Calculate volume averaged value of shape function derivative
  _avg_grad_test.resize(_test.size());
  for (_i = 0; _i < _test.size(); ++_i)
  {
    _avg_grad_test[_i] = 0.0;
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      _avg_grad_test[_i] += _grad_test[_i][_qp](_component) * _ad_JxW[_qp] * _ad_coord[_qp];

    _avg_grad_test[_i] /= ad_current_elem_volume;
  }
}