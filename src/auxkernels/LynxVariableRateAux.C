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
