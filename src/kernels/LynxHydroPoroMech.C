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
