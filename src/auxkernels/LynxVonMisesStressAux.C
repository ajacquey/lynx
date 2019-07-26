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

#include "LynxVonMisesStressAux.h"

registerMooseObject("LynxApp", LynxVonMisesStressAux);

template <>
InputParameters
validParams<LynxVonMisesStressAux>()
{
  InputParameters params = validParams<LynxStressAuxBase>();
  params.addClassDescription("Calculates the Von Mises stress.");
  return params;
}

LynxVonMisesStressAux::LynxVonMisesStressAux(const InputParameters & parameters)
  : LynxStressAuxBase(parameters)
{
}

Real
LynxVonMisesStressAux::computeValue()
{
  RankTwoTensor stress_dev = _stress[_qp].deviatoric();
  return std::sqrt(3.0 / 2.0) * stress_dev.L2norm();
}
