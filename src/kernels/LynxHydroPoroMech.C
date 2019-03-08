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

#include "LynxHydroPoroMech.h"
#include "MooseMesh.h"

registerMooseObject("LynxApp", LynxHydroPoroMech);

template <>
InputParameters
validParams<LynxHydroPoroMech>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Poromechanics coupling kernel.");
  params.addRequiredCoupledVar(
      "displacements", "The string of displacements variables suitable for the problem statement.");
  return params;
}

LynxHydroPoroMech::LynxHydroPoroMech(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _ndisp(coupledComponents("displacements")),
    _poro_mech(getMaterialProperty<Real>("poro_mechanical")),
    _poro_mech_jac(getMaterialProperty<Real>("poro_mechanical_jac")),
    _disp_var(_ndisp)
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError("The number of displacement variables supplied must match the mesh dimension.");

  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

Real
LynxHydroPoroMech::computeQpResidual()
{
  return _poro_mech[_qp] * _test[_i][_qp];
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

Real
LynxHydroPoroMech::computeQpJacobian()
{
  return 0.0;
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

Real
LynxHydroPoroMech::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    if (jvar == _disp_var[i])
      return _poro_mech_jac[_qp] * _grad_phi[_j][_qp](i) * _test[_i][_qp] / _dt;

  return 0.0;
}
