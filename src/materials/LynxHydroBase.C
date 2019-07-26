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

#include "LynxHydroBase.h"

template <>
InputParameters
validParams<LynxHydroBase>()
{
  InputParameters params = validParams<LynxMaterialBase>();
  params.addClassDescription("Base class for calculating the thermal properties.");
  params.addRequiredCoupledVar("porosity", "The porosity auxiliary variable.");
  return params;
}

LynxHydroBase::LynxHydroBase(const InputParameters & parameters)
  : LynxMaterialBase(parameters),
    _porosity(coupledValue("porosity")),
    _K(getDefaultMaterialProperty<Real>("bulk_modulus")),
    _tangent_modulus(getDefaultMaterialProperty<RankFourTensor>("tangent_modulus")),
    _strain_increment(getDefaultMaterialProperty<RankTwoTensor>("strain_increment")),
    _viscous_strain_incr(getDefaultMaterialProperty<RankTwoTensor>("viscous_strain_increment")),
    _plastic_strain_incr(getDefaultMaterialProperty<RankTwoTensor>("plastic_strain_increment")),
    _biot(declareProperty<Real>("biot_coefficient")),
    _C_d(declareProperty<Real>("bulk_compressibility")),
    _C_biot(declareProperty<Real>("biot_compressibility")),
    _fluid_mobility(declareProperty<Real>("fluid_mobility")),
    _poro_mech(declareProperty<Real>("poro_mechanical")),
    _poro_mech_jac(declareProperty<Real>("poro_mechanical_jac")),
    _C_f(_fe_problem.getMaxQps()),
    _C_s(_fe_problem.getMaxQps()),
    _k(_fe_problem.getMaxQps()),
    _eta_f(_fe_problem.getMaxQps())
{
}

void
LynxHydroBase::computeQpProperties()
{
  computeQpCompressibilities();
  computeQpFluidMobility();
  computeQpPoroMech();
}

void
LynxHydroBase::computeQpCompressibilities()
{
  computeQpFluidCompressibility();
  computeQpSolidCompressibility();

  // Drained compressibility
  _C_d[_qp] = (_K[_qp] != 0.0) ? 1.0 / _K[_qp] : 0.0;

  // Biot coefficient
  _biot[_qp] = 1.0;
  if (_C_d[_qp] != 0.0)
    _biot[_qp] -= _C_s[_qp] / _C_d[_qp];

  // Pore compressibility
  Real C_phi = (_biot[_qp] - _porosity[_qp]) * _C_d[_qp];

  // Biot compressibility
  _C_biot[_qp] = _porosity[_qp] * _C_f[_qp] + (1.0 - _biot[_qp]) * C_phi;
}

void
LynxHydroBase::computeQpFluidMobility()
{
  computeQpPermeability();
  computeQpFluidViscosity();

  // Fluid mobility
  _fluid_mobility[_qp] = _k[_qp] / _eta_f[_qp];
  if (_C_biot[_qp] != 0.0)
    _fluid_mobility[_qp] /= _C_biot[_qp];
}

void
LynxHydroBase::computeQpPoroMech()
{
  Real K_cto = _tangent_modulus[_qp].sum3x3() / 9.0;
  RankTwoTensor e_tot = _strain_increment[_qp] / _dt;
  RankTwoTensor e_in = (_viscous_strain_incr[_qp] + _plastic_strain_incr[_qp]) / _dt;
  _poro_mech[_qp] = _biot[_qp] * e_tot.trace() + (1.0 - _biot[_qp]) * e_in.trace();
  _poro_mech_jac[_qp] = _biot[_qp] + (1.0 - _biot[_qp]) * (1.0 - K_cto / _K[_qp]);

  if (_C_biot[_qp] != 0.0)
  {
    _poro_mech[_qp] /= _C_biot[_qp];
    _poro_mech_jac[_qp] /= _C_biot[_qp];
  }
}
