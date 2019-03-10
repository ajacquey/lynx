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

#include "LynxDensityConstant.h"

registerMooseObject("LynxApp", LynxDensityConstant);

template <>
InputParameters
validParams<LynxDensityConstant>()
{
  InputParameters params = validParams<LynxDensityBase>();
  params.addClassDescription("Material calculating densities as constant values.");
  return params;
}

LynxDensityConstant::LynxDensityConstant(const InputParameters & parameters)
  : LynxDensityBase(parameters)
{
}

void
LynxDensityConstant::computeQpProperties()
{
  computeQpGravity();
  _rho_f[_qp] = averageProperty(_fluid_density);
  _rho_s[_qp] = averageProperty(_solid_density);
  _rho_b[_qp] = _porosity[_qp] * _rho_f[_qp] + (1.0 - _porosity[_qp]) * _rho_s[_qp];
  _reference_rho_b[_qp] = _rho_b[_qp];
}
