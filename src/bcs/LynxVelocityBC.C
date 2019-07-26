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
