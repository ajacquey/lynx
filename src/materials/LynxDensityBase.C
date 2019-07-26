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

#include "LynxDensityBase.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<LynxDensityBase>()
{
  InputParameters params = validParams<LynxMaterialBase>();
  params.addClassDescription("Base class for calculating densities and gravity.");
  params.addCoupledVar("porosity", "The porosity auxiliary variable.");
  params.addParam<bool>("has_gravity", false, "Model with gravity on?");
  params.addParam<Real>("gravity_acceleration", 9.81, "The magnitude of the gravity acceleration.");
  params.addParam<std::vector<Real>>("fluid_density", "The fluid density.");
  params.addParam<std::vector<Real>>("solid_density", "The solid density.");
  return params;
}

LynxDensityBase::LynxDensityBase(const InputParameters & parameters)
  : LynxMaterialBase(parameters),
    _porosity(isCoupled("porosity") ? coupledValue("porosity") : _zero),
    _has_gravity(getParam<bool>("has_gravity")),
    _g(_has_gravity ? getParam<Real>("gravity_acceleration") : 0.0),
    _fluid_density(isParamValid("fluid_density") ? getLynxParam<Real>("fluid_density")
                                                 : std::vector<Real>(_n_composition, 0.0)),
    _solid_density(isParamValid("solid_density") ? getLynxParam<Real>("solid_density")
                                                 : std::vector<Real>(_n_composition, 0.0)),
    _gravity(declareProperty<RealVectorValue>("gravity_vector")),
    _rho_f(declareProperty<Real>("fluid_density")),
    _rho_s(declareProperty<Real>("solid_density")),
    _rho_b(declareProperty<Real>("bulk_density")),
    _reference_rho_b(declareProperty<Real>("reference_bulk_density"))
{
}

void
LynxDensityBase::computeQpGravity()
{
  if (_mesh.dimension() == 3)
    _gravity[_qp] = RealVectorValue(0., 0., -_g);
  else if (_mesh.dimension() == 2)
    _gravity[_qp] = RealVectorValue(0., -_g, 0.);
  else if (_mesh.dimension() == 1)
    _gravity[_qp] = RealVectorValue(-_g, 0., 0.);
}
