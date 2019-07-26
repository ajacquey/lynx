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
