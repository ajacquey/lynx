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

#include "LynxVelocityBC.h"
#include "Function.h"

registerMooseObject("LynxApp", LynxVelocityBC);

template <>
InputParameters
validParams<LynxVelocityBC>()
{
  InputParameters params = validParams<PresetNodalBC>();
  params.addClassDescription("Applies a velocity whose value is described by a function.");
  params.addParam<Real>("value", 0.0, "Value of the velocity applied.");
  params.addParam<FunctionName>("function", "Function giving the velocity applied.");
  return params;
}

LynxVelocityBC::LynxVelocityBC(const InputParameters & parameters)
  : PresetNodalBC(parameters),
    _u_old(valueOld()),
    _value(getParam<Real>("value")),
    _function(isParamValid("function") ? &getFunction("function") : NULL)
{
}

Real
LynxVelocityBC::computeQpValue()
{
  Real value = _value;

  if (_function)
    value = _function->value(_t, *_current_node);

  return _u_old[_qp] + value * _dt;
}