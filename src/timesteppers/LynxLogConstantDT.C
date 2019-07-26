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

#include "LynxLogConstantDT.h"

registerMooseObject("MooseApp", LynxLogConstantDT);

template <>
InputParameters
validParams<LynxLogConstantDT>()
{
  InputParameters params = validParams<TimeStepper>();
  params.addClassDescription(
      "TimeStepper which imposes a time step constant in the logarithmic space");
  params.addRequiredRangeCheckedParam<Real>("log_dt", "log_dt > 0", "Time step in log10(time)");
  params.addRequiredRangeCheckedParam<Real>(
      "first_dt", "first_dt > 0", "Initial time step (in absolute time)");
  params.addRequiredRangeCheckedParam<Real>("max_dt", "max_dt > 0", "Maximum value of time step.");
  params.addRangeCheckedParam<Real>(
      "growth_factor",
      2,
      "growth_factor>=1",
      "Maximum ratio of new to previous timestep sizes following a step that required the time"
      " step to be cut due to a failed solve.");
  return params;
}

LynxLogConstantDT::LynxLogConstantDT(const InputParameters & parameters)
  : TimeStepper(parameters),
    _log_dt(getParam<Real>("log_dt")),
    _first_dt(getParam<Real>("first_dt")),
    _max_dt(getParam<Real>("max_dt")),
    _dt_factor(std::pow(10.0, _log_dt)),
    _growth_factor(getParam<Real>("growth_factor"))
{
}

Real
LynxLogConstantDT::computeInitialDT()
{
  return _first_dt;
}

Real
LynxLogConstantDT::computeDT()
{
  Real next = _time * _dt_factor;
  return std::min(next - _time, std::min(_growth_factor * getCurrentDT(), _max_dt));
}
