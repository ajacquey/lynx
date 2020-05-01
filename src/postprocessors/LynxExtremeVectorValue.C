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

#include "LynxExtremeVectorValue.h"
#include "MooseMesh.h"

#include <algorithm>
#include <limits>

registerMooseObject("LynxApp", LynxExtremeVectorValue);

InputParameters
LynxExtremeVectorValue::validParams()
{
  InputParameters params = ElementExtremeValue::validParams();
  params.addClassDescription(
      "Compute the global minimum/maximum of a vectorial variable @ quadrature "
      "points with reference to the node.");
  params.addCoupledVar("add_var_1", 0.0, "The first additional variable.");
  params.addCoupledVar("add_var_2", 0.0, "The second additional variable.");
  return params;
}

LynxExtremeVectorValue::LynxExtremeVectorValue(const InputParameters & parameters)
  : ElementExtremeValue(parameters),
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
