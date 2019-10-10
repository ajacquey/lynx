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
#include "Function.h"

template <ComputeStage>
class LynxADDensityThermal;

declareADValidParams(LynxADDensityThermal);

template <ComputeStage compute_stage>
class LynxADDensityThermal : public LynxADDensityBase<compute_stage>
{
public:
  LynxADDensityThermal(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const ADVariableValue & _temp;
  const std::vector<Real> _beta_fluid;
  const std::vector<Real> _beta_solid;
  Real _temp_ref;
  const Function * _temp_ref_fct;

  usingDensityBaseMembers;
};