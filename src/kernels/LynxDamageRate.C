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

#include "LynxDamageRate.h"

registerMooseObject("LynxApp", LynxDamageRate);

template <>
InputParameters
validParams<LynxDamageRate>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Damage rate for damage rheology.");
  params.addCoupledVar("damage_rate", "The damage rate auxiliary variable");
  return params;
}

LynxDamageRate::LynxDamageRate(const InputParameters & parameters)
  : Kernel(parameters),
    _u_old(valueOld()),
    _coupled_dam(isCoupled("damage_rate")),
    _damage_rate(_coupled_dam ? coupledValue("damage_rate") : _zero)
{
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

Real
LynxDamageRate::computeQpResidual()
{
  Real rate = std::min(_damage_rate[_qp], (1.0 - _u_old[_qp]) / _dt);
  return -rate * _test[_i][_qp];
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

Real
LynxDamageRate::computeQpJacobian()
{
  return 0.0;
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

Real
LynxDamageRate::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
