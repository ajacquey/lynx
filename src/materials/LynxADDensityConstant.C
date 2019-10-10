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

#include "LynxADDensityConstant.h"

registerADMooseObject("LynxApp", LynxADDensityConstant);

defineADValidParams(
    LynxADDensityConstant,
    LynxADDensityBase,
  params.addClassDescription("Material calculating densities as constant values."););

template <ComputeStage compute_stage>
LynxADDensityConstant<compute_stage>::LynxADDensityConstant(const InputParameters & parameters)
  : LynxADDensityBase<compute_stage>(parameters)
{
}

template <ComputeStage compute_stage>
void
LynxADDensityConstant<compute_stage>::computeQpProperties()
{
  this->computeQpGravity();
  _rho_f[_qp] = this->averageProperty(_fluid_density);
  _rho_s[_qp] = this->averageProperty(_solid_density);
  if (_coupled_porosity)
    _rho_b[_qp] = (*_porosity)[_qp] * _rho_f[_qp] + (1.0 - (*_porosity)[_qp]) * _rho_s[_qp];
  else
    _rho_b[_qp] = _rho_s[_qp];
  _reference_rho_b[_qp] = _rho_b[_qp];
}

adBaseClass(LynxADDensityConstant);