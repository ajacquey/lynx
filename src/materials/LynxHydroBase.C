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
    _inelastic_strain(getDefaultMaterialProperty<RankTwoTensor>("inelastic_strain")),
    _inelastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("inelastic_strain")),
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
  RankTwoTensor e_in = (_inelastic_strain[_qp] - _inelastic_strain_old[_qp]) / _dt;
  _poro_mech[_qp] = _biot[_qp] * e_tot.trace() + (1.0 - _biot[_qp]) * e_in.trace();
  _poro_mech_jac[_qp] = _biot[_qp] + (1.0 - _biot[_qp]) * (1.0 - K_cto / _K[_qp]);

  if (_C_biot[_qp] != 0.0)
  {
    _poro_mech[_qp] /= _C_biot[_qp];
    _poro_mech_jac[_qp] /= _C_biot[_qp];
  }
}