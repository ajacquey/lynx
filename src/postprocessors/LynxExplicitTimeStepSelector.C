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

#include "LynxExplicitTimeStepSelector.h"
#include "libmesh/quadrature.h"

registerMooseObject("LynxApp", LynxExplicitTimeStepSelector);

template <>
InputParameters
validParams<LynxExplicitTimeStepSelector>()
{
  InputParameters params = validParams<ElementPostprocessor>();
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
  : DerivativeMaterialInterface<ElementPostprocessor>(parameters),
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
