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

#include "LynxADHydroPoroMech.h"

registerADMooseObject("LynxApp", LynxADHydroPoroMech);

defineADValidParams(LynxADHydroPoroMech,
                    ADKernel,
  params.addClassDescription("Poromechanics coupling kernel."););

template <ComputeStage compute_stage>
LynxADHydroPoroMech<compute_stage>::LynxADHydroPoroMech(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _poro_mech(getADMaterialProperty<Real>("poro_mechanical"))
{
}

template <ComputeStage compute_stage>
ADReal
LynxADHydroPoroMech<compute_stage>::computeQpResidual()
{
  return _poro_mech[_qp] * _test[_i][_qp];
}

adBaseClass(LynxADHydroPoroMech);