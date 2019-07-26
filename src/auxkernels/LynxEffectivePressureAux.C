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

#include "LynxEffectivePressureAux.h"

registerMooseObject("LynxApp", LynxEffectivePressureAux);

template <>
InputParameters
validParams<LynxEffectivePressureAux>()
{
  InputParameters params = validParams<LynxStressAuxBase>();
  params.addClassDescription("Calculates the effective pressure.");
  return params;
}

LynxEffectivePressureAux::LynxEffectivePressureAux(const InputParameters & parameters)
  : LynxStressAuxBase(parameters)
{
}

Real
LynxEffectivePressureAux::computeValue()
{
  return -_stress[_qp].trace() / 3.0;
}
