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

#include "LynxHeatFluxAux.h"

registerMooseObject("LynxApp", LynxHeatFluxAux);

template <>
InputParameters
validParams<LynxHeatFluxAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription(
    "Calculates the heat flux in each element for the given direction.");
  params.addRequiredCoupledVar("temperature", "The temperature variable.");
  params.addRequiredRangeCheckedParam<unsigned int>(
    "component",
    "component >= 0 & component <= 2",
    "Integer corresponding to direction of the heat flux (0, 1, 2.");
  return params;
}

LynxHeatFluxAux::LynxHeatFluxAux(const InputParameters & parameters)
  : DerivativeMaterialInterface<AuxKernel>(parameters),
    _component(getParam<unsigned int>("component")),
    _grad_T(coupledGradient("temperature")),
    _rhoC_b(getDefaultMaterialProperty<Real>("bulk_specific_heat")),
    _thermal_diff(getDefaultMaterialProperty<Real>("thermal_diffusivity"))
{
}

Real
LynxHeatFluxAux::computeValue()
{
  if (_rhoC_b[_qp] == 0.0)
    _cond_bulk = _thermal_diff[_qp];
  else
    _cond_bulk = _thermal_diff[_qp] * _rhoC_b[_qp];

  return -1*_cond_bulk*_grad_T[_qp](_component);
}