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

#include "LynxADAdvectionBase.h"
#include "Assembly.h"

InputParameters
LynxADAdvectionBase::validParams()
{
  InputParameters params = ADTimeKernel::validParams();
  params.addClassDescription("Base class for the advection term corrected by an artificial entropy "
                             "viscosity stabilization term after Guermond et al. (2011).");
  params.addRequiredCoupledVar("velocities", "The string of velocities that advect the problem.");
  params.addRequiredCoupledVar("entropy",
                               "The entropy used to calculate the artificial viscosity.");
  MooseEnum element_length_type("min=0 max=1 average=2", "min");
  params.addParam<MooseEnum>(
      "element_length_type", element_length_type, "The diameter of a single cell.");
  params.addParam<Real>("beta_stabilization", 0.026, "The beta local stabilization parameter.");
  params.addParam<Real>("cr_stabilization", 0.5, "The cr local stabilization parameter.");
  params.addParam<PostprocessorName>(
      "pp_max_vel", "The postprocessor to retrieve the maximum velocity on the whole domain.");
  params.addParam<PostprocessorName>(
      "pp_max_var",
      "The postprocessor to retrieve the maximum advected variable on the whole domain.");
  params.addParam<PostprocessorName>(
      "pp_min_var",
      "The postprocessor to retrieve the minimum advected variable on the whole domain.");
  params.addParam<PostprocessorName>(
      "pp_avg_var",
      "The postprocessor to retrieve the average advected variable on the whole domain.");
  params.addParam<PostprocessorName>(
      "pp_max_entropy", "The postprocessor to retrieve the maximum entropy on the whole domain.");
  params.addParam<PostprocessorName>(
      "pp_min_entropy", "The postprocessor to retrieve the minimum entropy on the whole domain.");
  params.addParam<PostprocessorName>(
      "pp_avg_entropy", "The postprocessor to retrieve the average entropy on the whole domain.");
  return params;
}

LynxADAdvectionBase::LynxADAdvectionBase(const InputParameters & parameters)
  : ADTimeKernel(parameters),
    _element_length_type(getParam<MooseEnum>("element_length_type")),
    _beta_stabilization(getParam<Real>("beta_stabilization")),
    _cr_stabilization(getParam<Real>("cr_stabilization")),
    _nvel(coupledComponents("velocities")),
    _vel(3),
    _vel_old(3),
    _vel_older(3),
    _value_old(this->_is_transient ? this->valueOld() : _zero),
    _gradient_old(this->_is_transient ? this->gradientOld() : _grad_zero),
    _second_old(this->_is_transient ? this->secondOld() : this->_second_zero),
    _value_older(this->_is_transient ? this->valueOlder() : _zero),
    _gradient_older(this->_is_transient ? this->gradientOlder() : _grad_zero),
    _second_older(this->_is_transient ? this->secondOlder() : this->_second_zero),
    _entropy_old(this->_is_transient ? this->coupledValueOld("entropy") : _zero),
    _entropy_older(this->_is_transient ? this->coupledValueOlder("entropy") : _zero),
    _pp_max_vel(getPostprocessorValue("pp_max_vel")),
    _pp_max_var(getPostprocessorValue("pp_max_var")),
    _pp_min_var(getPostprocessorValue("pp_min_var")),
    _pp_avg_var(getPostprocessorValue("pp_avg_var")),
    _pp_max_entropy(getPostprocessorValue("pp_max_entropy")),
    _pp_min_entropy(getPostprocessorValue("pp_min_entropy")),
    _pp_avg_entropy(getPostprocessorValue("pp_avg_entropy")),
    _residual(this->_fe_problem.getMaxQps())
{
  if (!this->_is_transient)
    mooseError("LynxADAdvectionBase: it does not make any sense to run this kernel on a "
               "steady state.");

  for (unsigned i = 0; i < _nvel; ++i)
  {
    _vel[i] = &adCoupledValue("velocities", i);
    _vel_old[i] = &coupledValueOld("velocities", i);
    _vel_older[i] = &coupledValueOlder("velocities", i);
  }
  for (unsigned i = _nvel; i < 3; ++i)
  {
    _vel[i] = &adZeroValue();
    _vel_old[i] = &_zero;
    _vel_older[i] = &_zero;
  }
}

Real
LynxADAdvectionBase::computeElementDiameter()
{
  Real diameter = 0.0;

  if (_current_elem->dim() == 1)
    diameter += _current_elem->volume();
  else
  {
    switch (_element_length_type)
    {
      case 0: // MIN
        diameter += _current_elem->hmin();
        break;
      case 1: // MAX
        diameter += _current_elem->hmax();
        break;
      case 2: // AVERAGE
        diameter += 0.5 * (_current_elem->hmin() + _current_elem->hmax());
        break;
      default:
        mooseError("LynxADAdvectionBase: unknown 'element_length_type'!");
    }
  }
  return diameter;
}

void
LynxADAdvectionBase::precalculateResidual()
{
  _artificial_viscosity = computeArtificialViscosity();
}

ADReal
LynxADAdvectionBase::computeQpResidual()
{
  ADRealVectorValue u((*_vel[0])[_qp], (*_vel[1])[_qp], (*_vel[2])[_qp]);
  ADReal res = u * _grad_u[_qp] * _test[_i][_qp];
  Real a2 = this->_t_step > 1 ? _dt / this->_dt_old : 0.0;
  Real a1 = 1.0 + a2;
  res += _artificial_viscosity *
         (_grad_test[_i][_qp] * (a1 * _gradient_old[_qp] - a2 * _gradient_older[_qp]));

  return res;
}