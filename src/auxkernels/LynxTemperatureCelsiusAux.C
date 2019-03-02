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

#include "LynxTemperatureCelsiusAux.h"

registerMooseObject("LynxApp", LynxTemperatureCelsiusAux);

template <>
InputParameters
validParams<LynxTemperatureCelsiusAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("temperature",
                               "The temperature variable to convert to Celsius.");
  params.addClassDescription("Converts Kelvin to degree Celsius.");
  return params;
}

LynxTemperatureCelsiusAux::LynxTemperatureCelsiusAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _T_K(coupledValue("temperature"))
{
}

Real
LynxTemperatureCelsiusAux::computeValue()
{
  Real dt_K = 273.15;
  return _T_K[_qp] - dt_K;
}
