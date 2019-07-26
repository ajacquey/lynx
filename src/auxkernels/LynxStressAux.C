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

#include "LynxStressAux.h"

registerMooseObject("LynxApp", LynxStressAux);

template <>
InputParameters
validParams<LynxStressAux>()
{
  InputParameters params = validParams<LynxStressAuxBase>();
  params.addClassDescription("Access a component of the stress tensor.");
  params.addCoupledVar("total_pressure", "The total pressure variable."); // For Stoke flow
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_i",
      "index_i >= 0 & index_i <= 2",
      "The index i of ij for the stress tensor (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_j",
      "index_j >= 0 & index_j <= 2",
      "The index j of ij for the stress tensor (0, 1, 2)");
  return params;
}

LynxStressAux::LynxStressAux(const InputParameters & parameters)
  : LynxStressAuxBase(parameters),
    _coupled_pt(isCoupled("total_pressure")),
    _pt(_coupled_pt ? coupledValue("total_pressure") : _zero),
    _i(getParam<unsigned int>("index_i")),
    _j(getParam<unsigned int>("index_j"))
{
}

Real
LynxStressAux::computeValue()
{
  return _stress[_qp](_i, _j) - (_i == _j) * _pt[_qp];
}
