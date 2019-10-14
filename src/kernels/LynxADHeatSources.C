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

#include "LynxADHeatSources.h"

registerADMooseObject("LynxApp", LynxADHeatSources);

defineADValidParams(
    LynxADHeatSources,
    ADKernel,
    params.addClassDescription(
        "Heat generation by radiogenic, shear heating and adiabatic processes.");
    params.addParam<Real>("coeff_shear_heating",
                          0.0,
                          "The coefficient in front of the shear heating generation.");
    params.addCoupledVar(
        "inelastic_heat",
        "The auxiliary variable holding the inelastic heat value for running in a subApp."););

template <ComputeStage compute_stage>
LynxADHeatSources<compute_stage>::LynxADHeatSources(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _coeff_Hs(getParam<Real>("coeff_shear_heating")),
    _rhoC_b(getADMaterialProperty<Real>("bulk_specific_heat")),
    _has_inelastic_heat_mat(hasMaterialProperty<Real>("inelastic_heat")),
    _radiogenic_heat(_has_inelastic_heat_mat ? &getADMaterialProperty<Real>("radiogenic_heat_production") : nullptr),
    _inelastic_heat_mat(_has_inelastic_heat_mat ? &getADMaterialProperty<Real>("inelastic_heat") : nullptr),
    _coupled_inelastic_heat(isCoupled("inelastic_heat")),
    _inelastic_heat(_coupled_inelastic_heat ? adCoupledValue("inelastic_heat") : adZeroValue())
{
}

template <ComputeStage compute_stage>
ADReal
LynxADHeatSources<compute_stage>::computeQpResidual()
{
  ADReal Hr = _has_inelastic_heat_mat ? (*_radiogenic_heat)[_qp] : 0.0;
  ADReal Hs = _coeff_Hs;
  if (_coupled_inelastic_heat)
    Hs *= _inelastic_heat[_qp];
  else if (_has_inelastic_heat_mat)
    Hs *= (*_inelastic_heat_mat)[_qp];
  else
    Hs *= 0.0;

  ADReal heat_sources = Hr + Hs;
  if (_rhoC_b[_qp] != 0.0)
    heat_sources /= _rhoC_b[_qp];

  return -heat_sources * _test[_i][_qp];
}