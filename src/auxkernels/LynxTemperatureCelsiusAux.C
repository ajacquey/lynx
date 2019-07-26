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

#include "LynxTemperatureCelsiusAux.h"

registerMooseObject("LynxApp", LynxTemperatureCelsiusAux);

template <>
InputParameters
validParams<LynxTemperatureCelsiusAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("temperature", "The temperature variable to convert to Celsius.");
  params.addClassDescription("Converts Kelvin to degree Celsius.");
  return params;
}

LynxTemperatureCelsiusAux::LynxTemperatureCelsiusAux(const InputParameters & parameters)
  : AuxKernel(parameters), _T_K(coupledValue("temperature"))
{
}

Real
LynxTemperatureCelsiusAux::computeValue()
{
  Real dt_K = 273.15;
  return _T_K[_qp] - dt_K;
}
