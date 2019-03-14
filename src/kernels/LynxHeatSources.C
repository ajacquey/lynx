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

#include "LynxHeatSources.h"

registerMooseObject("LynxApp", LynxHeatSources);

template <>
InputParameters
validParams<LynxHeatSources>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(
      "Heat generation by radiogenic, shear heating and adiabatic processes.");
  params.addParam<bool>(
      "use_displaced_mesh", true, "Set the displaced mesh flag to true by default.");
  params.addParam<Real>(
      "coeff_shear_heating", 0.0, "The coefficient in front of the shear heating generation.");
  params.addCoupledVar("displacements",
                       "The string of displacements suitable for the problem statement");
  params.addCoupledVar("inelastic_heat", "The auxiliary variable holding the inelastic heat value for running in a subApp.");
  return params;
}

LynxHeatSources::LynxHeatSources(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _coeff_Hs(getParam<Real>("coeff_shear_heating")),
    _rhoC_b(getDefaultMaterialProperty<Real>("bulk_specific_heat")),
    _radiogenic_heat(getDefaultMaterialProperty<Real>("radiogenic_heat_production")),
    _rho_b(getDefaultMaterialProperty<Real>("bulk_density")),
    _dinvrho_dtemp(getDefaultMaterialProperty<Real>("dinvrho_dtemp")),
    _inelastic_heat_mat(getDefaultMaterialProperty<Real>("inelastic_heat")),
    _adiabatic_heat(getDefaultMaterialProperty<Real>("adiabatic_heat")),
    _damage_heat(getDefaultMaterialProperty<Real>("damage_heat")),
    _dinelastic_heat_dstrain(getDefaultMaterialProperty<RankTwoTensor>("dinelastic_heat_dstrain")),
    _coupled_disp(isCoupled("displacements")),
    _ndisp(_coupled_disp ? coupledComponents("displacements") : 0),
    _disp_var(_ndisp),
    _coupled_inelastic_heat(isCoupled("inelastic_heat")),
    _inelastic_heat(_coupled_inelastic_heat ? coupledValue("inelastic_heat") : _zero)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);

  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_coupled_disp && (_ndisp != _mesh.dimension()))
    mooseError("LynxHeatSources: The number of displacement variables supplied must match the "
               "mesh dimension.");
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

Real
LynxHeatSources::computeQpResidual()
{
  Real Hr = _radiogenic_heat[_qp];
  Real Hs = _coeff_Hs;
  if (_coupled_inelastic_heat)
    Hs *= _inelastic_heat[_qp];
  else
    Hs *= _inelastic_heat_mat[_qp];
  Real Ha = _adiabatic_heat[_qp];
  Real Hd = _damage_heat[_qp];
  Real heat_sources = Hr + Hs + Ha + Hd;
  if (_rhoC_b[_qp] != 0.0)
    heat_sources /= _rhoC_b[_qp];

  return -heat_sources * _test[_i][_qp];
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

Real
LynxHeatSources::computeQpJacobian()
{
  Real Hr = _radiogenic_heat[_qp];
  Real Hs = _coeff_Hs;
  if (_coupled_inelastic_heat)
    Hs *= _inelastic_heat[_qp];
  else
    Hs *= _inelastic_heat_mat[_qp];
  Real Ha = _adiabatic_heat[_qp];
  Real Hd = _damage_heat[_qp];
  Real heat_sources = Hr + Hs + Ha + Hd;
  if (_rhoC_b[_qp] != 0.0)
    heat_sources /= _rhoC_b[_qp];

  return -heat_sources * _rho_b[_qp] * _dinvrho_dtemp[_qp] * _phi[_j][_qp] * _test[_i][_qp];
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

Real
LynxHeatSources::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    if (jvar == _disp_var[i])
      return -_test[_i][_qp] * _coeff_Hs * (_dinelastic_heat_dstrain[_qp] * _grad_phi[_j][_qp])(i);

  return 0.0;
}
