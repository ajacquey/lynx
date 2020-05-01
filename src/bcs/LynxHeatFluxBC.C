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

#include "LynxHeatFluxBC.h"
#include "Function.h"
#include "MooseRandom.h"

registerMooseObject("LynxApp", LynxHeatFluxBC);

InputParameters
LynxHeatFluxBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
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
