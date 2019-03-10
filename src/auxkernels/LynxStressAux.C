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
