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

#include "LynxElementAverageValue.h"

registerMooseObject("LynxApp", LynxElementAverageValue);

InputParameters
LynxElementAverageValue::validParams()
{
  InputParameters params = ElementAverageValue::validParams();
  params.addClassDescription("Compute the average value (integral sense) based on a forward "
                             "projection (unconditionally stable) of the advected variable.");
  return params;
}

LynxElementAverageValue::LynxElementAverageValue(const InputParameters & parameters)
  : ElementAverageValue(parameters),
    _value_old(_fe_problem.isTransient() ? valueOld() : _zero),
    _value_older(_fe_problem.isTransient() ? valueOlder() : _zero)
{
}

Real
LynxElementAverageValue::computeQpIntegral()
{
  return _t_step > 0 ? (1.0 + _dt / _dt_old) * _value_old[_qp] - _dt / _dt_old * _value_older[_qp]
                     : _value_old[_qp];
}
