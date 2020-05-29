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

#include "LynxADHydroDarcy.h"

registerMooseObject("LynxApp", LynxADHydroDarcy);

InputParameters
LynxADHydroDarcy::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Divergence of Darcy velocity kernel.");
  return params;
}

LynxADHydroDarcy::LynxADHydroDarcy(const InputParameters & parameters)
  : ADKernel(parameters),
    _fluid_mobility(getMaterialProperty<Real>("fluid_mobility")),
    _coupled_grav(hasMaterialProperty<RealVectorValue>("gravity_vector")),
    _gravity(_coupled_grav ? &getMaterialProperty<RealVectorValue>("gravity_vector") : nullptr),
    _rho_f(_coupled_grav ? &getADMaterialProperty<Real>("fluid_density") : nullptr)
{
}

ADReal
LynxADHydroDarcy::computeQpResidual()
{
  ADRealVectorValue grav_term = _coupled_grav ? -(*_rho_f)[_qp] * (*_gravity)[_qp] : ADRealVectorValue();

  return _fluid_mobility[_qp] * (_grad_u[_qp] + grav_term) * _grad_test[_i][_qp];
}