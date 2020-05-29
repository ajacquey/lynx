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

#include "LynxADDarcyVelocityAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LynxApp", LynxADDarcyVelocityAux);

InputParameters
LynxADDarcyVelocityAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the given component of the Darcy velocity.");
  params.addRequiredCoupledVar("fluid_pressure", "The fluid pressure variable.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "component",
      "component >= 0 & component <= 2",
      "Integer corresponding to direction of the Darcy velocity (0, 1, 2.");
  return params;
}

LynxADDarcyVelocityAux::LynxADDarcyVelocityAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _component(getParam<unsigned int>("component")),
    _grad_pf(coupledGradient("fluid_pressure")),
    _fluid_mobility(getMaterialProperty<Real>("fluid_mobility")),
    _coupled_grav(hasMaterialProperty<RealVectorValue>("gravity_vector")),
    _gravity(_coupled_grav ? &getMaterialProperty<RealVectorValue>("gravity_vector") : nullptr),
    _rho_f(_coupled_grav ? &getADMaterialProperty<Real>("fluid_density") : nullptr)
{
}

Real
LynxADDarcyVelocityAux::computeValue()
{
  RealVectorValue grav_term = _coupled_grav ? MetaPhysicL::raw_value(-(*_rho_f)[_qp]) *
                                                  (*_gravity)[_qp]
                                            : RealVectorValue();

  return -_fluid_mobility[_qp] *
         (_grad_pf[_qp](_component) + grav_term(_component));
}
