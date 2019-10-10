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

#include "LynxADAdvectionComposition.h"

registerADMooseObject("LynxApp", LynxADAdvectionComposition);

defineADValidParams(LynxADAdvectionComposition, LynxADAdvectionBase,
  params.addClassDescription("Composition advection kernel."););

template <ComputeStage compute_stage>
LynxADAdvectionComposition<compute_stage>::LynxADAdvectionComposition(const InputParameters & parameters)
  : LynxADAdvectionBase<compute_stage>(parameters)
{
}

template <ComputeStage compute_stage>
ADReal
LynxADAdvectionComposition<compute_stage>::computeArtificialViscosity()
{
  Real diameter = LynxADAdvectionBase<compute_stage>::computeElementDiameter();

  computeEntropyResidual();

  ADReal max_residual = 0.0;
  ADReal max_velocity = 0.0;

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    Real u0 = 0.5 * ((*_vel_old[0])[_qp] + (*_vel_older[0])[_qp]);
    Real u1 = 0.5 * ((*_vel_old[1])[_qp] + (*_vel_older[1])[_qp]);
    Real u2 = 0.5 * ((*_vel_old[2])[_qp] + (*_vel_older[2])[_qp]);

    RealVectorValue u(u0, u1, u2);
    max_residual = std::max(_residual[_qp], max_residual);
    max_velocity = std::max(u.norm(), max_velocity);
  }

  ADReal max_viscosity = _beta_stabilization * max_velocity * diameter;
  ADReal entropy_variation =
      std::max(_pp_max_entropy - _pp_avg_entropy, _pp_avg_entropy - _pp_min_entropy);

  if (this->_t_step <= 1 || std::abs(entropy_variation) < 1e-50)
    return max_viscosity;

  // If velocity is null, assume a sensible value to get back an artificial diffusion
  // here we do: velocity ~ _diffusivity / length_scale
  if (std::abs(_pp_max_vel) < 1e-50)
    return (_beta_stabilization * diameter) / this->_mesh.dimension();

  ADReal entropy_viscosity =
      _cr_stabilization * diameter * diameter * max_residual / entropy_variation;

  return std::min(max_viscosity, entropy_viscosity);
}

template <ComputeStage compute_stage>
void
LynxADAdvectionComposition<compute_stage>::computeEntropyResidual()
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

    _residual[_qp] = std::abs(dvar_dt + std::abs(field - _pp_avg_var) * u_grad_var);
  }
}
