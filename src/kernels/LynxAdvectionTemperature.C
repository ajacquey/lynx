/******************************************************************************/
/*                       LYNX, a MOOSE-based application                      */
/*                                                                            */
/*          Copyright (C) 2017 by Antoine B. Jacquey and Mauro Cacace         */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*    This program is free software: you can redistribute it and/or modify    */
/*    it under the terms of the GNU General Public License as published by    */
/*      the Free Software Foundation, either version 3 of the License, or     */
/*                     (at your option) any later version.                    */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program. If not, see <http://www.gnu.org/licenses/>     */
/******************************************************************************/

#include "LynxAdvectionTemperature.h"

registerMooseObject("LynxApp", LynxAdvectionTemperature);

template <>
InputParameters
validParams<LynxAdvectionTemperature>()
{
  InputParameters params = validParams<LynxAdvectionBase>();
  params.addClassDescription("Temperature advection kernel.");
  params.addParam<Real>(
      "coeff_shear_heating", 0.0, "The coefficient in front of the shear heating generation.");
  return params;
}

LynxAdvectionTemperature::LynxAdvectionTemperature(const InputParameters & parameters)
  : LynxAdvectionBase(parameters),
    _thermal_diff(getDefaultMaterialProperty<Real>("thermal_diffusivity")),
    _coeff_Hs(getParam<Real>("coeff_shear_heating")),
    _rhoC(getDefaultMaterialProperty<Real>("bulk_specific_heat")),
    _radiogenic_heat(getDefaultMaterialProperty<Real>("radiogenic_heat_production")),
    _inelastic_heat(getDefaultMaterialProperty<Real>("inelastic_heat"))
{
}

Real
LynxAdvectionTemperature::computeArtificialViscosity()
{
  Real diameter = computeElementDiameter();

  computeEntropyResidual();

  Real max_residual = 0.0;
  Real max_velocity = 0.0;
  Real max_diffusivity = 0.0;
  Real max_rho_cp = 0.0;

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

  // If velocity is null, assume a sensible value to get back an artificial diffusion
  // here we do: velocity ~ _diffusivity / length_scale
  if (std::abs(_pp_max_vel) < 1e-50)
    return _beta_stabilization * (max_diffusivity / max_rho_cp) / _mesh.dimension() * diameter;

  Real max_viscosity = _beta_stabilization * max_rho_cp * max_velocity * diameter;
  Real entropy_variation =
      std::max(_pp_max_entropy - _pp_avg_entropy, _pp_avg_entropy - _pp_min_entropy);

  if (_t_step <= 1 || std::abs(entropy_variation) < 1e-50)
    return max_viscosity;

  Real entropy_viscosity =
      _cr_stabilization * diameter * diameter * max_residual / entropy_variation;

  return std::min(max_viscosity, entropy_viscosity);
}

void
LynxAdvectionTemperature::computeEntropyResidual()
{
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    Real u0 = 0.5 * ((*_vel_old[0])[_qp] + (*_vel_older[0])[_qp]);
    Real u1 = 0.5 * ((*_vel_old[1])[_qp] + (*_vel_older[1])[_qp]);
    Real u2 = 0.5 * ((*_vel_old[2])[_qp] + (*_vel_older[2])[_qp]);
    RealVectorValue u(u0, u1, u2);

    Real dvar_dt = (_dt_old == 0.0) ? 0.0 : (_entropy_old[_qp] - _entropy_older[_qp]) / _dt_old;
    Real u_grad_var = u * (_gradient_old[_qp] + _gradient_older[_qp]) * 0.5;

    Real field = (_value_old[_qp] + _value_older[_qp]) * 0.5;

    Real laplace = 0.5 * (_second_old[_qp].tr() + _second_older[_qp].tr());
    Real kappa_laplace_var = _thermal_diff[_qp] * laplace;

    Real Hr = _radiogenic_heat[_qp];
    Real Hs = _coeff_Hs * _inelastic_heat[_qp];
    Real heat_sources = Hr + Hs;
    if (_rhoC[_qp] != 0.0)
      heat_sources /= _rhoC[_qp];

    _residual[_qp] = std::abs(dvar_dt + std::abs(field - _pp_avg_var) * u_grad_var -
                              kappa_laplace_var - heat_sources);
  }
}
