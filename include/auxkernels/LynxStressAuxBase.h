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

#pragma once

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "DerivativeMaterialInterface.h"

class LynxStressAuxBase;

template <>
InputParameters validParams<LynxStressAuxBase>();

class LynxStressAuxBase : public DerivativeMaterialInterface<AuxKernel>
{
public:
  LynxStressAuxBase(const InputParameters & parameters);

protected:
  const MaterialProperty<RankTwoTensor> & _stress;
};