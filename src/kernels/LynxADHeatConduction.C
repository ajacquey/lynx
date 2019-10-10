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

#include "LynxADHeatConduction.h"

registerADMooseObject("LynxApp", LynxADHeatConduction);

defineADValidParams(
    LynxADHeatConduction,
    ADKernel,
  params.addClassDescription("Heat conduction kernel.");
  params.addParam<bool>(
      "use_displaced_mesh", true, "Set the displaced mesh flag to true by default."););

template <ComputeStage compute_stage>
LynxADHeatConduction<compute_stage>::LynxADHeatConduction(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _thermal_diff(getADMaterialProperty<Real>("thermal_diffusivity"))
{
}

template <ComputeStage compute_stage>
ADReal
LynxADHeatConduction<compute_stage>::computeQpResidual()
{
  return _thermal_diff[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
}

adBaseClass(LynxADHeatConduction);