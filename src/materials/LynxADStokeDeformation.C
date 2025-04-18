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

#include "LynxADStokeDeformation.h"

registerMooseObject("LynxApp", LynxADStokeDeformation);

InputParameters
LynxADStokeDeformation::validParams()
{
  InputParameters params = LynxADDeformationBase::validParams();
  params.addClassDescription("Class calculating strain and stress for a viscous (Stoke) rheology.");
  params.addCoupledVar("dynamic_pressure", "The dynamic pressure variable.");
  // Stoke parameters
  params.addParam<std::vector<Real>>("A_diffusion",
                                     "List of pre-exponential factors for diffusion creep.");
  params.addParam<std::vector<Real>>("E_diffusion",
                                     "List of activation energy for diffusion creep.");
  params.addParam<std::vector<Real>>("V_diffusion",
                                     "List of activation (molar) volume for diffusion creep.");
  params.addParam<std::vector<Real>>("A_dislocation",
                                     "List of pre-exponential factors for dislocation creep.");
  params.addParam<std::vector<Real>>("E_dislocation",
                                     "List of activation energy for dislocation creep.");
  params.addParam<std::vector<Real>>("V_dislocation",
                                     "List of activation (molar) volume for dislocation creep.");
  params.addParam<std::vector<Real>>("n_dislocation",
                                     "List of power law exponents for dislocation creep.");
  params.addParam<Real>("gas_constant", 8.3144621, "The universal gas constant");
  params.addParam<std::vector<Real>>(
      "initial_viscosity", "The vector of initial viscosity for purely viscous deformation.");
  params.addParam<Real>("background_strain_rate",
                        "The background strain rate for purely viscous deformation.");
  params.addParam<std::vector<Real>>(
      "eta_min", "The lower threshold for the effective viscosity for purely viscous deformation.");
  params.addParam<std::vector<Real>>(
      "eta_max", "The upper threshold for the effective viscosity for purely viscous deformation.");
  return params;
}

LynxADStokeDeformation::LynxADStokeDeformation(const InputParameters & parameters)
  : LynxADDeformationBase(parameters),
    _pdyn(isCoupled("dynamic_pressure") ? adCoupledValue("dynamic_pressure") : _ad_zero),
    _temp(_coupled_temp ? adCoupledValue("temperature") : _ad_zero),
    // Stoke parameters
    _has_diffusion_creep(isParamValid("A_diffusion")),
    _A_diffusion(_has_diffusion_creep ? getLynxParam<Real>("A_diffusion")
                                      : std::vector<Real>(_n_composition, 0.0)),
    _E_diffusion((_has_diffusion_creep && isParamValid("E_diffusion"))
                     ? getLynxParam<Real>("E_diffusion")
                     : std::vector<Real>(_n_composition, 0.0)),
    _V_diffusion((_has_diffusion_creep && isParamValid("V_diffusion"))
                     ? getLynxParam<Real>("V_diffusion")
                     : std::vector<Real>(_n_composition, 0.0)),
    _has_dislocation_creep(isParamValid("A_dislocation")),
    _A_dislocation(_has_dislocation_creep ? getLynxParam<Real>("A_dislocation")
                                          : std::vector<Real>(_n_composition, 0.0)),
    _n_dislocation((_has_dislocation_creep && isParamValid("n_dislocation"))
                       ? getLynxParam<Real>("n_dislocation")
                       : std::vector<Real>(_n_composition, 1.0)),
    _E_dislocation((_has_dislocation_creep && isParamValid("E_dislocation"))
                       ? getLynxParam<Real>("E_dislocation")
                       : std::vector<Real>(_n_composition, 0.0)),
    _V_dislocation((_has_dislocation_creep && isParamValid("V_dislocation"))
                       ? getLynxParam<Real>("V_dislocation")
                       : std::vector<Real>(_n_composition, 0.0)),
    _gas_constant(getParam<Real>("gas_constant")),
    _has_background_strain_rate(isParamValid("background_strain_rate")),
    _has_initial_viscosity(_has_background_strain_rate ? false
                                                       : isParamValid("_has_initial_viscosity")),
    _initial_viscosity(_has_initial_viscosity ? getLynxParam<Real>("initial_viscosity")
                                              : std::vector<Real>(_n_composition, 0.0)),
    _background_strain_rate(_has_background_strain_rate ? getParam<Real>("background_strain_rate")
                                                        : 0.0),
    _eta_min(isParamValid("eta_min") ? getLynxParam<Real>("eta_min")
                                     : std::vector<Real>(_n_composition, 0.0)),
    _eta_max(isParamValid("eta_max") ? getLynxParam<Real>("eta_max")
                                     : std::vector<Real>(_n_composition, 1.0e+99)),
    // Stoke properties
    _eta_eff(declareADProperty<Real>("effective_viscosity"))
{
  // Check consistency
  if (!_has_diffusion_creep && !_has_dislocation_creep)
    mooseError("LynxStokeDeformation: you need to provide either a diffusion or a dislocation "
               "creep model!");

  // Check if initial viscosity
  if (_has_dislocation_creep && !_has_diffusion_creep && !_has_initial_viscosity &&
      !_has_background_strain_rate)
    mooseWarning("LynxStokeDeformation: you provided a dislocation creep model but no initial "
                 "viscosity or background strain rate!");
}

void
LynxADStokeDeformation::initializeQpDeformation()
{
  _A_diff = averageProperty(_A_diffusion);
  _E_diff = averageProperty(_E_diffusion);
  _V_diff = averageProperty(_V_diffusion);

  _A_disl = averageProperty(_A_dislocation);
  _n_disl = averageProperty(_n_dislocation);
  _E_disl = averageProperty(_E_dislocation);
  _V_disl = averageProperty(_V_dislocation);
}

ADReal
LynxADStokeDeformation::volumetricDeformation()
{
  return _pdyn[_qp] + _plith[_qp];
}

ADRankTwoTensor
LynxADStokeDeformation::deviatoricDeformation(const ADReal & pressure)
{
  // Compute effective viscosity
  _eta_eff[_qp] = computeQpEffectiveViscosity(pressure);

  return 2.0 * _eta_eff[_qp] * _strain_increment[_qp].deviatoric() / _dt;
}

ADReal
LynxADStokeDeformation::computeQpEffectiveViscosity(const ADReal & pressure)
{
  ADReal strain_rate_II = std::sqrt(0.5) * _strain_increment[_qp].deviatoric().L2norm() / _dt;

  if (_t_step <= 1 && strain_rate_II == 0.0 && _has_initial_viscosity)
    return std::min(std::max(averageProperty(_initial_viscosity), averageProperty(_eta_min)),
                    averageProperty(_eta_max));

  strain_rate_II = (_has_background_strain_rate && _t_step <= 1 && strain_rate_II == 0.0)
                       ? _background_strain_rate
                       : strain_rate_II;

  ADReal one_on_eta_diff = 0.0, one_on_eta_disl = 0.0;
  if (_coupled_temp)
  {
    ADReal RT = _gas_constant * _temp[_qp];
    _A_diff *= std::exp(-(_E_diff + pressure * _V_diff) / RT);
    _A_disl *= std::exp(-(_E_disl + pressure * _V_disl) / RT);
  }

  if (_has_diffusion_creep)
    one_on_eta_diff = computeQpOneOnDiffViscosity(_A_diff);
  if (_has_dislocation_creep)
    one_on_eta_disl = computeQpOneOnDislViscosity(_A_disl, _n_disl, strain_rate_II);

  ADReal eta = 1.0 / (one_on_eta_diff + one_on_eta_disl);

  return std::min(std::max(eta, averageProperty(_eta_min)), averageProperty(_eta_max));
}

ADReal
LynxADStokeDeformation::computeQpOneOnDiffViscosity(const ADReal A)
{
  return 2.0 * A;
}

ADReal
LynxADStokeDeformation::computeQpOneOnDislViscosity(const ADReal A,
                                                    const ADReal n,
                                                    const ADReal eII)
{
  if ((eII == 0.0) && (n == 1.0))
    return 2.0;
  else
    return 2.0 * std::pow(A, 1.0 / n) * std::pow(eII, 1.0 - 1.0 / n);
}

void
LynxADStokeDeformation::computeQpThermalSources()
{
  _inelastic_heat[_qp] =
      _stress[_qp].deviatoric().doubleContraction(_strain_increment[_qp].deviatoric() / _dt);
}