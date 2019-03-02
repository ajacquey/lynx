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

#include "LynxFunctionNoiseIC.h"
#include "Function.h"
#include "MooseRandom.h"

registerMooseObject("LynxApp", LynxFunctionNoiseIC);

template <>
InputParameters
validParams<LynxFunctionNoiseIC>()
{
  InputParameters params = validParams<InitialCondition>();
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