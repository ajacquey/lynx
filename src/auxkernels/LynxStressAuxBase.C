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

#include "LynxStressAuxBase.h"

template <>
InputParameters
validParams<LynxStressAuxBase>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Base class for outputting stress values.");
  return params;
}

LynxStressAuxBase::LynxStressAuxBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<AuxKernel>(parameters),
    _stress(getMaterialProperty<RankTwoTensor>("stress"))
{
}
