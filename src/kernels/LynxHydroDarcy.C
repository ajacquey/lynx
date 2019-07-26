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
