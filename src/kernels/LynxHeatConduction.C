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
