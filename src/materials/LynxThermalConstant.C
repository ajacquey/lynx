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

#include "LynxThermalConstant.h"

registerMooseObject("LynxApp", LynxThermalConstant);

template <>
InputParameters
validParams<LynxThermalConstant>()
{
  InputParameters params = validParams<LynxThermalBase>();
  params.addClassDescription("Constant thermal properties.");
  params.addParam<std::vector<Real>>("fluid_thermal_conductivity",
                                     "The fluid thermal conductivity.");
  params.addRequiredParam<std::vector<Real>>("solid_thermal_conductivity",
                                             "The solid thermal conductivity.");
  params.addParam<std::vector<Real>>("fluid_heat_capacity", "The fluid heat capacity.");
  params.addParam<std::vector<Real>>("solid_heat_capacity", "The solid heat capacity.");
  params.addParam<std::vector<Real>>("fluid_thermal_expansion", "The fluid volumetric thermal expansion coefficient.");
  params.addParam<std::vector<Real>>("solid_thermal_expansion", "The solid volumetric thermal expansion coefficient.");
  return params;
}

LynxThermalConstant::LynxThermalConstant(const InputParameters & parameters)
  : LynxThermalBase(parameters),
    _fluid_thermal_cond(isParamValid("fluid_thermal_conductivity")
                            ? getLynxParam<Real>("fluid_thermal_conductivity")
                            : std::vector<Real>(_n_composition, 0.0)),
    _solid_thermal_cond(getLynxParam<Real>("solid_thermal_conductivity")),
    _fluid_heat_cap(isParamValid("fluid_heat_capacity") ? getLynxParam<Real>("fluid_heat_capacity")
                                                        : std::vector<Real>(_n_composition, 0.0)),
    _solid_heat_cap(isParamValid("solid_heat_capacity") ? getLynxParam<Real>("solid_heat_capacity")
                                                        : std::vector<Real>(_n_composition, 0.0)),
    _fluid_thermal_exp(isParamValid("fluid_thermal_expansion")
                           ? getLynxParam<Real>("fluid_thermal_expansion")
                           : std::vector<Real>(_n_composition, 0.0)),
    _solid_thermal_exp(isParamValid("solid_thermal_expansion")
                           ? getLynxParam<Real>("solid_thermal_expansion")
                           : std::vector<Real>(_n_composition, 0.0))
{
}

void
LynxThermalConstant::computeQpHeatCap()
{
  _c_f[_qp] = averageProperty(_fluid_heat_cap);
  _c_s[_qp] = averageProperty(_solid_heat_cap);
}

void
LynxThermalConstant::computeQpThermalCond()
{
  _lambda_f[_qp] = averageProperty(_fluid_thermal_cond);
  _lambda_s[_qp] = averageProperty(_solid_thermal_cond);
}

void
LynxThermalConstant::computeQpThermalExp()
{
  _beta_f[_qp] = averageProperty(_fluid_thermal_exp);
  _beta_s[_qp] = averageProperty(_solid_thermal_exp);
}