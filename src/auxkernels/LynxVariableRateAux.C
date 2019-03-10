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

#include "LynxVariableRateAux.h"

registerMooseObject("LynxApp", LynxVariableRateAux);

template <>
InputParameters
validParams<LynxVariableRateAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Computes the rate of change of a coupled scalar variable.");
  params.addRequiredCoupledVar("coupled_variable", "The coupled variable.");
  params.addParam<Real>("time_scale_factor", 1.0, "Multiply the time step size by this factor.");
  params.addParam<bool>(
      "output_relative",
      false,
      "If true, write gradient relative to the absolute value of coupled_variable.");
  return params;
}

LynxVariableRateAux::LynxVariableRateAux(const InputParameters & parameters)
  : DerivativeMaterialInterface<AuxKernel>(parameters),
    _cvar(coupledValue("coupled_variable")),
    _cvar_old(coupledValueOld("coupled_variable")),
    _tscale(getParam<Real>("time_scale_factor")),
    _relative(getParam<bool>("output_relative"))
{
}

Real
LynxVariableRateAux::computeValue()
{
  if (_relative)
    return (_cvar[_qp] - _cvar_old[_qp]) / _cvar_old[_qp] / _dt / _tscale;
  else
    return (_cvar[_qp] - _cvar_old[_qp]) / _dt / _tscale;
}
