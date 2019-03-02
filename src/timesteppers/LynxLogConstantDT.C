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
