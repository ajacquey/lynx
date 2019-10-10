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
class LynxADPressureLoad;

declareADValidParams(LynxADPressureLoad);

template <ComputeStage compute_stage>
class LynxADPressureLoad : public ADKernel<compute_stage>
{
public:
  LynxADPressureLoad(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty(Real) & _bulk_density;
  const ADMaterialProperty(RealVectorValue) & _gravity;

  usingKernelMembers;
};