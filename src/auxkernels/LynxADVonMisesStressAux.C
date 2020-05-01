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

#include "LynxADVonMisesStressAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LynxApp", LynxADVonMisesStressAux);

InputParameters
LynxADVonMisesStressAux::validParams()
{
  InputParameters params = LynxADStressAuxBase::validParams();
  params.addClassDescription("Calculates the Von Mises stress.");
  return params;
}

LynxADVonMisesStressAux::LynxADVonMisesStressAux(const InputParameters & parameters)
  : LynxADStressAuxBase(parameters)
{
}

Real
LynxADVonMisesStressAux::computeValue()
{
  ADRankTwoTensor stress_dev = _stress[_qp].deviatoric();
  return std::sqrt(3.0 / 2.0) * MetaPhysicL::raw_value(stress_dev.L2norm());
}
