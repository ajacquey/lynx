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

#include "LynxADDensityBase.h"

template <ComputeStage>
class LynxADDensityConstant;

declareADValidParams(LynxADDensityConstant);

template <ComputeStage compute_stage>
class LynxADDensityConstant : public LynxADDensityBase<compute_stage>
{
public:
  LynxADDensityConstant(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  usingDensityBaseMembers;
};