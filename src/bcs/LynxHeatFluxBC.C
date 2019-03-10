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

#include "LynxHeatFluxBC.h"
#include "Function.h"
#include "MooseRandom.h"

registerMooseObject("LynxApp", LynxHeatFluxBC);

template <>
InputParameters
validParams<LynxHeatFluxBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addClassDescription("Applies a constant heat flux.");
  params.addParam<Real>("value", 0.0, "Value of the heat flux applied.");
  params.addParam<FunctionName>("function", "Function giving the heat flux applied.");
  params.addRangeCheckedParam<Real>("random_percentage",
                                    0.0,
                                    "random_percentage>=0 & random_percentage<=1",
                                    "The percentage (0 to 1) defining the random range.");
  return params;
}

LynxHeatFluxBC::LynxHeatFluxBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _value(getParam<Real>("value")),
    _function(isParamValid("function") ? &getFunction("function") : NULL),
    _rand_per(getParam<Real>("random_percentage")),
    _rhoC_b(getMaterialProperty<Real>("bulk_specific_heat"))
{
}

Real
LynxHeatFluxBC::computeQpResidual()
{
  Real value = _value;

  if (_function)
    value = _function->value(_t, _q_point[_qp]);

  Real min = (1.0 - _rand_per) * value;
  Real range = 2.0 * _rand_per * value;

  Real rand_num = MooseRandom::rand();
  rand_num *= range;
  rand_num += min;

  if (_rhoC_b[_qp] != 0.0)
    rand_num /= _rhoC_b[_qp];

  return -_test[_i][_qp] * rand_num;
}
