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

#include "LynxDensityThermal.h"

registerMooseObject("LynxApp", LynxDensityThermal);

template <>
InputParameters
validParams<LynxDensityThermal>()
{
  InputParameters params = validParams<LynxDensityBase>();
  params.addClassDescription(
      "Material calculating densities as a simple linear function of temperature.");
  params.addParam<std::vector<Real>>("beta_fluid", "The fluid thermal expansion coefficient.");
  params.addParam<std::vector<Real>>("beta_solid", "The solid thermal expansion coefficient.");
  params.addParam<Real>("reference_temperature", 0.0, "The reference temperature.");
  params.addParam<FunctionName>("reference_temperature_fct",
                                "The reference temperature given by a function.");
  return params;
}

LynxDensityThermal::LynxDensityThermal(const InputParameters & parameters)
  : LynxDensityBase(parameters),
    _temperature(coupledValue("temperature")),
    _drho_dtemp(declareProperty<Real>("drho_dtemp")),
    _dinvrho_dtemp(declareProperty<Real>("dinvrho_dtemp")),
    _beta_fluid(isParamValid("beta_fluid") ? getLynxParam<Real>("beta_fluid")
                                           : std::vector<Real>(_n_composition, 0.0)),
    _beta_solid(isParamValid("beta_solid") ? getLynxParam<Real>("beta_solid")
                                           : std::vector<Real>(_n_composition, 0.0)),
    _temp_ref(getParam<Real>("reference_temperature")),
    _temp_ref_fct(isParamValid("reference_temperature_fct")
                      ? &getFunction("reference_temperature_fct")
                      : NULL)
{
}

void
LynxDensityThermal::computeQpProperties()
{
  computeQpGravity();

  Real temp_ref = _temp_ref;
  if (_temp_ref_fct)
    temp_ref = _temp_ref_fct->value(_t, _q_point[_qp]);

  _rho_f[_qp] = averageProperty(_fluid_density) *
                (1.0 - averageProperty(_beta_fluid) * (_temperature[_qp] - temp_ref));
  _rho_s[_qp] = averageProperty(_solid_density) *
                (1.0 - averageProperty(_beta_solid) * (_temperature[_qp] - temp_ref));

  _rho_b[_qp] = _porosity[_qp] * _rho_f[_qp] + (1.0 - _porosity[_qp]) * _rho_s[_qp];
  _reference_rho_b[_qp] = _rho_b[_qp];

  Real drho_f_dtemp = -1.0 * averageProperty(_beta_fluid) * averageProperty(_fluid_density);
  Real drho_s_dtemp = -1.0 * averageProperty(_beta_solid) * averageProperty(_solid_density);

  _drho_dtemp[_qp] = _porosity[_qp] * drho_f_dtemp + (1.0 - _porosity[_qp]) * drho_s_dtemp;

  _dinvrho_dtemp[_qp] =
      (_rho_b[_qp] > 0.0) ? -1.0 / std::pow(_rho_b[_qp], 2.0) * _drho_dtemp[_qp] : 0.0;
}
