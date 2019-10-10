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

#include "LynxADDensityBase.h"
#include "MooseMesh.h"

defineADValidParams(LynxADDensityBase,
                    LynxADMaterialBase,
  params.addClassDescription("Base class for calculating densities and gravity.");
  params.addCoupledVar("porosity","The porosity auxiliary variable.");
  params.addParam<bool>("has_gravity", false, "Model with gravity on?");
  params.addParam<Real>("gravity_acceleration", 9.81, "The magnitude of the gravity acceleration.");
  params.addParam<std::vector<Real>>("fluid_density", "The fluid density.");
  params.addParam<std::vector<Real>>("solid_density", "The solid density."););

template <ComputeStage compute_stage>
LynxADDensityBase<compute_stage>::LynxADDensityBase(const InputParameters & parameters)
  : LynxADMaterialBase<compute_stage>(parameters),
    _porosity(isCoupled("porosity") ? adCoupledValue("porosity") : adZeroValue()),
    _has_gravity(getParam<bool>("has_gravity")),
    _g(_has_gravity ? getParam<Real>("gravity_acceleration") : 0.0),
    _fluid_density(isParamValid("fluid_density") ? this->getLynxParam("fluid_density")
                                                 : std::vector<Real>(_n_composition, 0.0)),
    _solid_density(isParamValid("solid_density") ? this->getLynxParam("solid_density")
                                                 : std::vector<Real>(_n_composition, 0.0)),
    _gravity(declareADProperty<RealVectorValue>("gravity_vector")),
    _rho_f(declareADProperty<Real>("fluid_density")),
    _rho_s(declareADProperty<Real>("solid_density")),
    _rho_b(declareADProperty<Real>("bulk_density")),
    _reference_rho_b(declareADProperty<Real>("reference_bulk_density"))
{
}

template <ComputeStage compute_stage>
void
LynxADDensityBase<compute_stage>::computeQpGravity()
{
  if (_mesh.dimension() == 3)
    _gravity[_qp] = ADRealVectorValue(0., 0., -_g);
  else if (_mesh.dimension() == 2)
    _gravity[_qp] = ADRealVectorValue(0., -_g, 0.);
  else if (_mesh.dimension() == 1)
    _gravity[_qp] = ADRealVectorValue(-_g, 0., 0.);
}

adBaseClass(LynxADDensityBase);