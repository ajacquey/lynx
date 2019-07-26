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

#include "LynxAdvectionBase.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<LynxAdvectionBase>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Base class for the advection term corrected by an artificial entropy "
                             "viscosity stabilization term after Guermond et al. (2011).");
  params.addCoupledVar("displacements", "The string of displacements that advect the problem.");
  params.addRequiredCoupledVar("velocities", "The string of velocities that advect the problem.");
  params.addCoupledVar("entropy", "The entropy used to calculate the artificial viscosity.");
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

LynxAdvectionBase::LynxAdvectionBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _element_length_type(getParam<MooseEnum>("element_length_type")),
    _beta_stabilization(getParam<Real>("beta_stabilization")),
    _cr_stabilization(getParam<Real>("cr_stabilization")),
    _has_disp(isCoupled("displacements")),
    _nvel(coupledComponents("velocities")),
    _disp_var(_nvel),
    _vel(3),
    _vel_old(3),
    _vel_older(3),
    _value_old(_fe_problem.isTransient() ? valueOld() : _zero),
    _gradient_old(_fe_problem.isTransient() ? gradientOld() : _grad_zero),
    _second_old(_fe_problem.isTransient() ? secondOld() : _second_zero),
    _value_older(_fe_problem.isTransient() ? valueOlder() : _zero),
    _gradient_older(_fe_problem.isTransient() ? gradientOlder() : _grad_zero),
    _second_older(_fe_problem.isTransient() ? secondOlder() : _second_zero),
    _entropy_old(_fe_problem.isTransient() ? coupledValueOld("entropy") : _zero),
    _entropy_older(_fe_problem.isTransient() ? coupledValueOlder("entropy") : _zero),
    _pp_max_vel(getPostprocessorValue("pp_max_vel")),
    _pp_max_var(getPostprocessorValue("pp_max_var")),
    _pp_min_var(getPostprocessorValue("pp_min_var")),
    _pp_avg_var(getPostprocessorValue("pp_avg_var")),
    _pp_max_entropy(getPostprocessorValue("pp_max_entropy")),
    _pp_min_entropy(getPostprocessorValue("pp_min_entropy")),
    _pp_avg_entropy(getPostprocessorValue("pp_avg_entropy")),
    _residual(_fe_problem.getMaxQps())
{
  if (!_fe_problem.isTransient())
    mooseError("LynxAdvectionBase: it does not make any sense to run this kernel on a "
               "steady state.");
  if ((_nvel != _mesh.dimension()))
    mooseError("LynxAdvectionBase: the number of supplied 'displacements' does not match the "
               "problem dimension.");
  if (_has_disp && (_nvel != coupledComponents("displacements")))
    mooseError("LynxAdvectionBase: the number of supplied 'displacements' should be the same as "
               "the number of supplied 'velocities'!");

  for (unsigned i = 0; i < _nvel; ++i)
  {
    _disp_var[i] = _has_disp ? coupled("displacements", i) : -1;
    _vel[i] = &coupledValue("velocities", i);
    _vel_old[i] = &coupledValueOld("velocities", i);
    _vel_older[i] = &coupledValueOlder("velocities", i);
  }
  for (unsigned i = _nvel; i < 3; ++i)
  {
    _vel[i] = &_zero;
    _vel_old[i] = &_zero;
    _vel_older[i] = &_zero;
  }
}

Real
LynxAdvectionBase::computeElementDiameter()
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
        mooseError("LynxAdvectionBase: unknow 'element_length_type'!");
    }
  }
  return diameter;
}

/******************************************************************************
 *                                   RESIDUALS                                 *
 *******************************************************************************/
void
LynxAdvectionBase::precalculateResidual()
{
  _artificial_viscosity = computeArtificialViscosity();
}

Real
LynxAdvectionBase::computeQpResidual()
{
  RealVectorValue u((*_vel[0])[_qp], (*_vel[1])[_qp], (*_vel[2])[_qp]);
  Real res = u * _grad_u[_qp] * _test[_i][_qp];
  Real a2 = _t_step > 1 ? _dt / _dt_old : 0.0;
  Real a1 = 1.0 + a2;
  res += _artificial_viscosity *
         (_grad_test[_i][_qp] * (a1 * _gradient_old[_qp] - a2 * _gradient_older[_qp]));

  return res;
}

/******************************************************************************
 *                                  JACOBIAN                                   *
 *******************************************************************************/
Real
LynxAdvectionBase::computeQpJacobian()
{
  RealVectorValue u((*_vel[0])[_qp], (*_vel[1])[_qp], (*_vel[2])[_qp]);
  Real jac = u * _grad_phi[_j][_qp] * _test[_i][_qp];

  return jac;
}

/******************************************************************************
 *                               OFF-DIAG JACOBIAN                             *
 *******************************************************************************/
Real
LynxAdvectionBase::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real jac = 0.0;
  for (unsigned int coupled_component = 0; coupled_component < _nvel; ++coupled_component)
    if (_has_disp && (jvar == _disp_var[coupled_component]))
      jac += _phi[_j][_qp] / _dt * _grad_u[_qp](coupled_component) * _test[_i][_qp];

  return jac;
}
