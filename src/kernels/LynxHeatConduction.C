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

#include "LynxHeatConduction.h"

registerMooseObject("LynxApp", LynxHeatConduction);

template <>
InputParameters
validParams<LynxHeatConduction>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Heat conduction kernel.");
  params.addParam<bool>(
      "use_displaced_mesh", true, "Set the displaced mesh flag to true by default.");
  return params;
}

LynxHeatConduction::LynxHeatConduction(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _thermal_diff(getMaterialProperty<Real>("thermal_diffusivity")),
    _rho(getDefaultMaterialProperty<Real>("bulk_density")),
    _dinvrho_dtemp(getDefaultMaterialProperty<Real>("dinvrho_dtemp"))
{
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

Real
LynxHeatConduction::computeQpResidual()
{
  return _thermal_diff[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

Real
LynxHeatConduction::computeQpJacobian()
{
  Real jac = _thermal_diff[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
  jac += (_thermal_diff[_qp] * _rho[_qp]) * _dinvrho_dtemp[_qp] * _phi[_j][_qp] * _grad_u[_qp] *
         _grad_test[_i][_qp];

  return jac;
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

// Real
// LynxHeatConduction::computeQpOffDiagJacobian(unsigned int jvar)
// {
//   return 0.0;
// }
