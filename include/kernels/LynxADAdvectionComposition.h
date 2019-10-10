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

#include "LynxADAdvectionBase.h"

template <ComputeStage>
class LynxADAdvectionComposition;

declareADValidParams(LynxADAdvectionComposition);

template <ComputeStage compute_stage>
class LynxADAdvectionComposition : public LynxADAdvectionBase<compute_stage>
{
public:
  LynxADAdvectionComposition(const InputParameters & parameters);

protected:
  virtual ADReal computeArtificialViscosity() override;
  virtual void computeEntropyResidual();

  usingAdvectionBaseMembers;
};