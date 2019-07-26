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

#include "LynxExtrapolateValue.h"

#include <algorithm>
#include <limits>

registerMooseObject("LynxApp", LynxExtrapolateValue);

template <>
InputParameters
validParams<LynxExtrapolateValue>()
{
  InputParameters params = validParams<ElementExtremeValue>();
  params.addClassDescription("Compute the global range of variation based on a forward "
                             "projection (unconditionally stable) of the advected variable.");
  return params;
}

LynxExtrapolateValue::LynxExtrapolateValue(const InputParameters & parameters)
  : DerivativeMaterialInterface<ElementExtremeValue>(parameters),
    _value_old(_fe_problem.isTransient() ? valueOld() : _zero),
    _value_older(_fe_problem.isTransient() ? valueOlder() : _zero)
{
}

void
LynxExtrapolateValue::computeQpValue()
{
  Real extrapolated_value =
      _t_step > 0 ? (1.0 + _dt / _dt_old) * _value_old[_qp] - _dt / _dt_old * _value_older[_qp]
                  : _value_old[_qp];
  switch (_type)
  {
    case MAX:
      _value = std::max(_value, extrapolated_value);
      break;
    case MIN:
      _value = std::min(_value, extrapolated_value);
      break;
  }
}
