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

#include "LynxADEffectivePressureAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LynxApp", LynxADEffectivePressureAux);

InputParameters
LynxADEffectivePressureAux::validParams()
{
  InputParameters params = LynxADStressAuxBase::validParams();
  params.addClassDescription("Calculates the effective pressure.");
  return params;
}

LynxADEffectivePressureAux::LynxADEffectivePressureAux(const InputParameters & parameters)
  : LynxADStressAuxBase(parameters)
{
}

Real
LynxADEffectivePressureAux::computeValue()
{
  return -MetaPhysicL::raw_value(_stress[_qp].trace()) / 3.0;
}
