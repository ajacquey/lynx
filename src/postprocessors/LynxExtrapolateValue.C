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
