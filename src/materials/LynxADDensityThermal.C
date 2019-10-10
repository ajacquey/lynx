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

#include "LynxADDensityThermal.h"

registerADMooseObject("LynxApp", LynxADDensityThermal);

defineADValidParams(
    LynxADDensityThermal,
    LynxADDensityBase,
  params.addClassDescription("Material calculating densities as a simple linear function of temperature.");
  params.addRequiredCoupledVar("temperature", "The temperature variable.");
  params.addParam<std::vector<Real>>("beta_fluid", "The fluid thermal expansion coefficient.");
  params.addParam<std::vector<Real>>("beta_solid", "The solid thermal expansion coefficient.");
  params.addParam<Real>("reference_temperature", 0.0, "The reference temperature.");
  params.addParam<FunctionName>("reference_temperature_fct",
                                "The reference temperature given by a function."););

template <ComputeStage compute_stage>
LynxADDensityThermal<compute_stage>::LynxADDensityThermal(const InputParameters & parameters)
  : LynxADDensityBase<compute_stage>(parameters),
    _temp(adCoupledValue("temperature")),
    _beta_fluid(isParamValid("beta_fluid") ? this->getLynxParam("beta_fluid")
                                           : std::vector<Real>(_n_composition, 0.0)),
    _beta_solid(isParamValid("beta_solid") ? this->getLynxParam("beta_solid")
                                           : std::vector<Real>(_n_composition, 0.0)),
    _temp_ref(getParam<Real>("reference_temperature")),
    _temp_ref_fct(isParamValid("reference_temperature_fct")
                      ? &getFunction("reference_temperature_fct")
                      : nullptr)
{
}

template <ComputeStage compute_stage>
void
LynxADDensityThermal<compute_stage>::computeQpProperties()
{
  this->computeQpGravity();

  Real temp_ref = _temp_ref;
  if (_temp_ref_fct)
    temp_ref = _temp_ref_fct->value(_t, _q_point[_qp]);

  _rho_f[_qp] = this->averageProperty(_fluid_density) *
                (1.0 - this->averageProperty(_beta_fluid) * (_temp[_qp] - temp_ref));
  _rho_s[_qp] = this->averageProperty(_solid_density) *
                (1.0 - this->averageProperty(_beta_solid) * (_temp[_qp] - temp_ref));

  if (_coupled_porosity)
    _rho_b[_qp] = (*_porosity)[_qp] * _rho_f[_qp] + (1.0 - (*_porosity)[_qp]) * _rho_s[_qp];
  else
    _rho_b[_qp] = _rho_s[_qp];

  _reference_rho_b[_qp] = _rho_b[_qp];
}

adBaseClass(LynxADDensityThermal);