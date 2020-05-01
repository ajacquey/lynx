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

#include "LynxExplicitTimeStepSelector.h"
#include "libmesh/quadrature.h"

registerMooseObject("LynxApp", LynxExplicitTimeStepSelector);

InputParameters
LynxExplicitTimeStepSelector::validParams()
{
  InputParameters params = ElementPostprocessor::validParams();
  params.addClassDescription("Postprocessor that computes the minimum value of h_min/|u|, where "
                             "|u| is coupled in as an aux variable.");
  params.addRequiredCoupledVar("vel_norm", "Velocity magnitude");
  params.addRequiredParam<Real>("beta",
                                "0 < beta < 1, choose some fraction of the limiting timestep size");
  params.addParam<Real>("epsilon", 0, "The epsilon parameter.");
  params.addParam<bool>(
      "has_premult",
      false,
      "Whether to premultiply the courant criterion by geometric and interpolation constraints.");
  params.addParam<Real>("initial_value", 1.0, "The initial time step size.");
  params.addParam<Real>("maximum_value", 1.0, "The maximum time step size.");
  return params;
}

LynxExplicitTimeStepSelector::LynxExplicitTimeStepSelector(const InputParameters & parameters)
  : ElementPostprocessor(parameters),
    _vel_norm(coupledValue("vel_norm")),
    _beta(getParam<Real>("beta")),
    _epsilon(isParamValid("epsilon") ? getParam<Real>("epsilon")
                                     : std::numeric_limits<Real>::epsilon()),
    _has_premult(getParam<bool>("has_premult")),
    _initial_value(getParam<Real>("initial_value")),
    _maximum_value(getParam<Real>("maximum_value"))
{
}

LynxExplicitTimeStepSelector::~LynxExplicitTimeStepSelector() {}

void
LynxExplicitTimeStepSelector::initialize()
{
  _value = std::numeric_limits<Real>::max();
}

void
LynxExplicitTimeStepSelector::execute()
{
  Real h_min = _current_elem->hmin();
  Real dim = static_cast<Real>(_current_elem->dim());
  for (unsigned qp = 0; qp < _qrule->n_points(); ++qp)
  {
    Real vel_norm = std::max(_vel_norm[qp], _epsilon);
    Real premult =
        _has_premult ? 1.0 / (1.7 * dim * std::sqrt(1.0 * dim)) / _qrule->n_points() : 1.0;
    Real courant_limit_dt = premult * h_min / vel_norm;

    _value = std::min(_value, _beta * courant_limit_dt);
  }
  if (_t_step <= 1)
    _value = _initial_value;
  if (_value > _maximum_value)
    _value = _maximum_value;
}

Real
LynxExplicitTimeStepSelector::getValue()
{
  _communicator.min(_value);
  return _value;
}

void
LynxExplicitTimeStepSelector::threadJoin(const UserObject & uo)
{
  const LynxExplicitTimeStepSelector & pps = dynamic_cast<const LynxExplicitTimeStepSelector &>(uo);
  _value = std::min(_value, pps._value);
}
