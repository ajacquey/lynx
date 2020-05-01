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

#include "LynxADStressAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LynxApp", LynxADStressAux);

InputParameters
LynxADStressAux::validParams()
{
  InputParameters params = LynxADStressAuxBase::validParams();
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

LynxADStressAux::LynxADStressAux(const InputParameters & parameters)
  : LynxADStressAuxBase(parameters),
    _coupled_pt(isCoupled("total_pressure")),
    _pt(_coupled_pt ? adCoupledValue("total_pressure") : _ad_zero),
    _i(getParam<unsigned int>("index_i")),
    _j(getParam<unsigned int>("index_j"))
{
}

Real
LynxADStressAux::computeValue()
{
  return MetaPhysicL::raw_value(_stress[_qp](_i, _j)) -
         (_i == _j) * MetaPhysicL::raw_value(_pt[_qp]);
}
