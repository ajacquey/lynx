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

#include "ADKernel.h"

template <ComputeStage>
class LynxADMass;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;

declareADValidParams(LynxADMass);

template <ComputeStage compute_stage>
class LynxADMass : public ADKernel<compute_stage>
{
public:
  LynxADMass(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _penalty;
  const unsigned int _penalty_type;
  const ADMaterialProperty(RankTwoTensor) & _strain_increment;

  usingKernelMembers;
};