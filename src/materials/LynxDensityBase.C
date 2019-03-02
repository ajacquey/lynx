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
  // params.addParam<std::vector<Real>>("fluid_density", std::vector<Real>(1, 0.0), "The fluid
  // density.");
  params.addParam<std::vector<Real>>("solid_density", "The solid density.");
  return params;
}

LynxDensityBase::LynxDensityBase(const InputParameters & parameters)
  : LynxMaterialBase(parameters),
    _porosity(isCoupled("porosity") ? coupledValue("porosity") : _zero),
    _has_gravity(getParam<bool>("has_gravity")),
    _g(_has_gravity ? getParam<Real>("gravity_acceleration") : 0.0),
    // _fluid_density(getLynxParam<Real>("fluid_density")),
    _solid_density(isParamValid("solid_density") ? getLynxParam<Real>("solid_density")
                                                 : std::vector<Real>(_n_composition, 0.0)),
    _gravity(declareProperty<RealVectorValue>("gravity_vector")),
    // _rho_f(declareProperty<Real>("fluid_density")),
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
