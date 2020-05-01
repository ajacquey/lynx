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

#include "LynxFunctionNoiseIC.h"
#include "Function.h"
#include "MooseRandom.h"

registerMooseObject("LynxApp", LynxFunctionNoiseIC);

InputParameters
LynxFunctionNoiseIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addClassDescription(
      "Applies a initial condition based on a function and adds some white noise.");
  params.addRequiredParam<FunctionName>("function", "The initial condition function.");
  params.addRequiredRangeCheckedParam<Real>("random_percentage",
                                            "random_percentage>0 & random_percentage<1",
                                            "The percentage (0 to 1) defining the random range.");
  return params;
}

LynxFunctionNoiseIC::LynxFunctionNoiseIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _func(getFunction("function")),
    _rand_per(getParam<Real>("random_percentage"))
{
}

Real
LynxFunctionNoiseIC::value(const Point & p)
{
  Real func_value = _func.value(_t, p);
  Real min = (1.0 - _rand_per) * func_value;
  Real range = 2.0 * _rand_per * func_value;

  Real rand_num = MooseRandom::rand();
  rand_num *= range;
  rand_num += min;

  return rand_num;
}

RealGradient
LynxFunctionNoiseIC::gradient(const Point & p)
{
  return _func.gradient(_t, p);
}