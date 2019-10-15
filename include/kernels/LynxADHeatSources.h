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
class LynxADHeatSources;

declareADValidParams(LynxADHeatSources);

template <ComputeStage compute_stage>
class LynxADHeatSources : public ADKernel<compute_stage>
{
public:
  LynxADHeatSources(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _coeff_Hs;
  unsigned int _nvel;
  std::vector<const ADVariableValue *> _vel;
  const ADVariableGradient & _grad_pressure; 
  const ADMaterialProperty(Real) & _rhoC_b;
  const bool _has_inelastic_heat_mat;
  const ADMaterialProperty(Real) * _radiogenic_heat;
  const ADMaterialProperty(Real) * _inelastic_heat_mat;
  const bool _coupled_inelastic_heat;
  const ADVariableValue & _inelastic_heat;
  const ADMaterialProperty(Real) & _thermal_exp;

  usingKernelMembers;
};