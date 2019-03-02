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

#include "LynxPressureLoad.h"

registerMooseObject("LynxApp", LynxPressureLoad);

template <>
InputParameters
validParams<LynxPressureLoad>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(
      "Kernel calculating the lithostatic pressure based on density distribution.");
  return params;
}

LynxPressureLoad::LynxPressureLoad(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _bulk_density(getDefaultMaterialProperty<Real>("reference_bulk_density")),
    _gravity(getDefaultMaterialProperty<RealVectorValue>("gravity_vector"))
{
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

Real
LynxPressureLoad::computeQpResidual()
{
  return (_grad_u[_qp] - _bulk_density[_qp] * _gravity[_qp]) * _grad_test[_i][_qp];
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

Real
LynxPressureLoad::computeQpJacobian()
{
  return _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
