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

#include "LynxDarcyVelocityAux.h"

registerMooseObject("LynxApp", LynxDarcyVelocityAux);

template <>
InputParameters
validParams<LynxDarcyVelocityAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculates the given component of the Darcy velocity.");
  params.addRequiredCoupledVar("fluid_pressure", "The fluid pressure variable.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "component",
      "component >= 0 & component <= 2",
      "Integer corresponding to direction of the Darcy velocity (0, 1, 2.");
  return params;
}

LynxDarcyVelocityAux::LynxDarcyVelocityAux(const InputParameters & parameters)
  : DerivativeMaterialInterface<AuxKernel>(parameters),
    _component(getParam<unsigned int>("component")),
    _grad_pf(coupledGradient("fluid_pressure")),
    _fluid_mobility(getDefaultMaterialProperty<Real>("fluid_mobility")),
    _rho_f(getDefaultMaterialProperty<Real>("fluid_density")),
    _gravity(getDefaultMaterialProperty<RealVectorValue>("gravity_vector"))    
{
}

Real
LynxDarcyVelocityAux::computeValue()
{
  return -_fluid_mobility[_qp] * (_grad_pf[_qp](_component) - _rho_f[_qp] * _gravity[_qp](_component));
}
