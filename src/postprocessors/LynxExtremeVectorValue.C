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

#include "LynxExtremeVectorValue.h"
#include "MooseMesh.h"

#include <algorithm>
#include <limits>

registerMooseObject("LynxApp", LynxExtremeVectorValue);

template <>
InputParameters
validParams<LynxExtremeVectorValue>()
{
  InputParameters params = validParams<ElementExtremeValue>();
  params.addClassDescription(
      "Compute the global minimum/maximum of a vectorial variable @ quadrature "
      "points with reference to the node.");
  params.addCoupledVar("add_var_1", 0.0, "The first additional variable.");
  params.addCoupledVar("add_var_2", 0.0, "The second additional variable.");
  return params;
}

LynxExtremeVectorValue::LynxExtremeVectorValue(const InputParameters & parameters)
  : DerivativeMaterialInterface<ElementExtremeValue>(parameters),
    _v(_mesh.dimension() > 1 ? coupledValue("add_var_1") : _zero),
    _w(_mesh.dimension() > 2 ? coupledValue("add_var_2") : _zero)
{
}

void
LynxExtremeVectorValue::computeQpValue()
{
  RealVectorValue vv(_u[_qp], _v[_qp], _w[_qp]);
  switch (_type)
  {
    case MAX:
      _value = std::max(_value, vv.norm());
      break;
    case MIN:
      _value = std::min(_value, vv.norm());
      break;
  }
}
