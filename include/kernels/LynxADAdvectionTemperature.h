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
class LynxADAdvectionTemperature;

declareADValidParams(LynxADAdvectionTemperature);

template <ComputeStage compute_stage>
class LynxADAdvectionTemperature : public LynxADAdvectionBase<compute_stage>
{
public:
  LynxADAdvectionTemperature(const InputParameters & parameters);

protected:
  virtual ADReal computeArtificialViscosity() override;
  virtual void computeEntropyResidual();

  const Real _coeff_Hs;

  const ADMaterialProperty(Real) & _thermal_diff;
  const ADMaterialProperty(Real) & _rhoC;
  const ADMaterialProperty(Real) & _radiogenic_heat;
  const ADMaterialProperty(Real) & _inelastic_heat;

  usingAdvectionBaseMembers;
};