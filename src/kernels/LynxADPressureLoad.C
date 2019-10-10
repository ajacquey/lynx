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

#include "LynxADPressureLoad.h"

registerADMooseObject("LynxApp", LynxADPressureLoad);

defineADValidParams(
    LynxADPressureLoad,
    ADKernel,
    params.addClassDescription(
        "Kernel calculating the lithostatic pressure based on density distribution."););

template <ComputeStage compute_stage>
LynxADPressureLoad<compute_stage>::LynxADPressureLoad(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _bulk_density(getADMaterialProperty<Real>("reference_bulk_density")),
    _gravity(getADMaterialProperty<RealVectorValue>("gravity_vector"))
{
}

template <ComputeStage compute_stage>
ADReal
LynxADPressureLoad<compute_stage>::computeQpResidual()
{
  return (_grad_u[_qp] - _bulk_density[_qp] * _gravity[_qp]) * _grad_test[_i][_qp];
}