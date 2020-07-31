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

#include "LynxADAdvectionTemperature.h"

registerMooseObject("LynxApp", LynxADAdvectionTemperature);

InputParameters
LynxADAdvectionTemperature::validParams()
{
  InputParameters params = LynxADAdvectionBase::validParams();
  params.addClassDescription("Temperature advection kernel.");
  params.addParam<Real>(
      "coeff_shear_heating", 0.0, "The coefficient in front of the shear heating generation.");
  params.addCoupledVar(
        "inelastic_heat",
        "The auxiliary variable holding the inelastic heat value for running in a subApp.");
  params.addCoupledVar("pressure", "The pressure variable");
  return params;
}

LynxADAdvectionTemperature::LynxADAdvectionTemperature(const InputParameters & parameters)
  : LynxADAdvectionBase(parameters),
    _coeff_Hs(getParam<Real>("coeff_shear_heating")),
    _grad_pressure(isCoupled("pressure") ? adCoupledGradient("pressure") : _ad_grad_zero),
    _thermal_diff(getADMaterialProperty<Real>("thermal_diffusivity")),
    _rhoC(getADMaterialProperty<Real>("bulk_specific_heat")),
    _has_inelastic_heat_mat(hasMaterialProperty<Real>("inelastic_heat")),
    _radiogenic_heat(_has_inelastic_heat_mat
                         ? &getADMaterialProperty<Real>("radiogenic_heat_production")
                         : nullptr),
    _inelastic_heat_mat(_has_inelastic_heat_mat ? &getADMaterialProperty<Real>("inelastic_heat")
                                                : nullptr),
    _coupled_inelastic_heat(isCoupled("inelastic_heat")),
    _inelastic_heat(_coupled_inelastic_heat ? adCoupledValue("inelastic_heat") : _ad_zero),
    _thermal_exp(getADMaterialProperty<Real>("thermal_expansion_coefficient"))
{
}

ADReal
LynxADAdvectionTemperature::computeArtificialViscosity()
{
  Real diameter = LynxADAdvectionBase::computeElementDiameter();

  computeEntropyResidual();

  ADReal max_residual = 0.0;
  ADReal max_velocity = 0.0;
  ADReal max_diffusivity = 0.0;
  ADReal max_rho_cp = 0.0;

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    Real u0 = 0.5 * ((*_vel_old[0])[_qp] + (*_vel_older[0])[_qp]);
    Real u1 = 0.5 * ((*_vel_old[1])[_qp] + (*_vel_older[1])[_qp]);
    Real u2 = 0.5 * ((*_vel_old[2])[_qp] + (*_vel_older[2])[_qp]);
    RealVectorValue u(u0, u1, u2);

    max_residual = std::max(_residual[_qp], max_residual);
    max_velocity = std::max(u.norm(), max_velocity);
    max_diffusivity = std::max(_thermal_diff[_qp], max_diffusivity);
    max_rho_cp = std::max(_rhoC[_qp], max_rho_cp);
  }

  ADReal max_viscosity = _beta_stabilization * max_rho_cp * max_velocity * diameter;
  ADReal entropy_variation =
      std::max(_pp_max_entropy - _pp_avg_entropy, _pp_avg_entropy - _pp_min_entropy);
      
  if (this->_t_step <= 1 || std::abs(entropy_variation) < 1e-50)
    return max_viscosity;

  // If velocity is null, assume a sensible value to get back an artificial diffusion
  // here we do: velocity ~ _diffusivity / length_scale
  if (std::abs(_pp_max_vel) < 1e-50)
    return _beta_stabilization * (max_diffusivity / max_rho_cp) / this->_mesh.dimension() * diameter;

  ADReal entropy_viscosity =
      _cr_stabilization * diameter * diameter * max_residual / entropy_variation;

  return std::min(max_viscosity, entropy_viscosity);
}

void
LynxADAdvectionTemperature::computeEntropyResidual()
{
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    Real u0 = 0.5 * ((*_vel_old[0])[_qp] + (*_vel_older[0])[_qp]);
    Real u1 = 0.5 * ((*_vel_old[1])[_qp] + (*_vel_older[1])[_qp]);
    Real u2 = 0.5 * ((*_vel_old[2])[_qp] + (*_vel_older[2])[_qp]);
    RealVectorValue u(u0, u1, u2);

    Real dvar_dt = (this->_dt_old == 0.0) ? 0.0 : (_entropy_old[_qp] - _entropy_older[_qp]) / this->_dt_old;
    Real u_grad_var = u * (_gradient_old[_qp] + _gradient_older[_qp]) * 0.5;

    Real field = (_value_old[_qp] + _value_older[_qp]) * 0.5;

    Real laplace = 0.5 * (_second_old[_qp].tr() + _second_older[_qp].tr());
    ADReal kappa_laplace_var = _thermal_diff[_qp] * laplace;

    ADReal Hr = _has_inelastic_heat_mat ? (*_radiogenic_heat)[_qp] : 0.0;
    ADReal Hs = _coeff_Hs;
    if (_coupled_inelastic_heat)
      Hs *= _inelastic_heat[_qp];
    else if (_has_inelastic_heat_mat)
      Hs *= (*_inelastic_heat_mat)[_qp];
    else
      Hs *= 0.0;

    ADRealVectorValue vel((*_vel[0])[_qp], (*_vel[1])[_qp], (*_vel[2])[_qp]);
    ADReal Ha = _thermal_exp[_qp] * _u[_qp] * vel * _grad_pressure[_qp];

    ADReal heat_sources = Hr + Hs + Ha;

    if (_rhoC[_qp] != 0.0)
      heat_sources /= _rhoC[_qp];

    _residual[_qp] = std::abs(dvar_dt + std::abs(field - _pp_avg_var) * u_grad_var -
                              kappa_laplace_var - heat_sources);
  }
}