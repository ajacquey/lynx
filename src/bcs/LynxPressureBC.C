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

#include "LynxPressureBC.h"
#include "Function.h"

registerMooseObject("LynxApp", LynxPressureBC);

template <>
InputParameters
validParams<LynxPressureBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addClassDescription("Applies a pressure on a given boundary in a given direction.");
  params.addRequiredParam<unsigned int>("component", "The component for the pressure.");
  params.addParam<Real>("value", 0.0, "Value of the pressure applied.");
  params.addParam<FunctionName>("function", "The function that describes the pressure.");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

LynxPressureBC::LynxPressureBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _component(getParam<unsigned int>("component")),
    _value(getParam<Real>("value")),
    _function(isParamValid("function") ? &getFunction("function") : NULL)
{
  if (_component > 2)
    mooseError("Invalid component given for ", name(), ": ", _component, ".\n");
}

Real
LynxPressureBC::computeQpResidual()
{
  Real value = _value;

  if (_function)
    value = _function->value(_t, _q_point[_qp]);

  return value * (_normals[_qp](_component) * _test[_i][_qp]);
}