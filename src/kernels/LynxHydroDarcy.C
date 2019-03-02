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

#include "LynxHydroDarcy.h"

registerMooseObject("LynxApp", LynxHydroDarcy);

template <>
InputParameters
validParams<LynxHydroDarcy>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Divergence of Darcy velocity kernel.");
  return params;
}

LynxHydroDarcy::LynxHydroDarcy(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _fluid_mobility(getMaterialProperty<Real>("fluid_mobility")),
    _gravity(getDefaultMaterialProperty<RealVectorValue>("gravity_vector")),
    _rho_f(getDefaultMaterialProperty<Real>("fluid_density"))
{
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

Real
LynxHydroDarcy::computeQpResidual()
{
  RealVectorValue grav_term = -_rho_f[_qp] * _gravity[_qp];

  return _fluid_mobility[_qp] * (_grad_u[_qp] + grav_term) * _grad_test[_i][_qp];
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

Real
LynxHydroDarcy::computeQpJacobian()
{
  return _fluid_mobility[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

// Real
// LynxHydroDarcy::computeQpOffDiagJacobian(unsigned int jvar)
// {
//   return 0.0;
// }