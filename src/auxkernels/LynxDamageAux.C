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

#include "LynxDamageAux.h"

registerMooseObject("LynxApp", LynxDamageAux);

template <>
InputParameters
validParams<LynxDamageAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Updating the damage auxiliary variable.");
  return params;
}

LynxDamageAux::LynxDamageAux(const InputParameters & parameters)
  : DerivativeMaterialInterface<AuxKernel>(parameters),
    _damage_rate(getDefaultMaterialProperty<Real>("damage_rate"))
{
}

Real
LynxDamageAux::computeValue()
{
  Real damage_rate = _damage_rate[_qp];
  if (damage_rate > (1.0 - _u_old[_qp]) / _dt) 
    damage_rate = (1.0 - _u_old[_qp]) / _dt;

  return _u_old[_qp] + damage_rate * _dt;
}