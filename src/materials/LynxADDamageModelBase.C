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

#include "LynxADDamageModelBase.h"
#include "MooseRandom.h"

InputParameters
LynxADDamageModelBase::validParams()
{
  InputParameters params = LynxADMaterialBase::validParams();
  params.addClassDescription("Base class for the damage rheology.");
  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");
  params.addCoupledVar("fluid_pressure", 0, "The fluid pressure variable.");
  params.addRangeCheckedParam<Real>("abs_tolerance",
                                    1.0e-10,
                                    "abs_tolerance > 0.0",
                                    "The absolute tolerance for the iterative update.");
  params.addRangeCheckedParam<Real>("rel_tolerance",
                                    1.0e-10,
                                    "rel_tolerance > 0.0",
                                    "The relative tolerance for the iterative update.");
  params.addRangeCheckedParam<unsigned int>(
      "max_iterations",
      500,
      "max_iterations >= 1",
      "The maximum number of iterations for the iterative update");
  params.addParam<std::vector<Real>>("initial_damage", "The initial damage value.");
  params.addParam<std::vector<Real>>("initial_damage_random", "The random initial damage window.");
  params.addParam<std::vector<Real>>("reference_fluid_pressure", "The reference fluid pressure.");
  return params;
}

LynxADDamageModelBase::LynxADDamageModelBase(const InputParameters & parameters)
  : LynxADMaterialBase(parameters),
    _pf(adCoupledValue("fluid_pressure")),
    _abs_tol(getParam<Real>("abs_tolerance")),
    _rel_tol(getParam<Real>("rel_tolerance")),
    _max_its(getParam<unsigned int>("max_iterations")),
    _damage0(isParamValid("initial_damage") ? getLynxParam<Real>("initial_damage")
                                            : std::vector<Real>(_n_composition, 0.0)),
    _damage0_rand(isParamValid("initial_damage_random") ? getLynxParam<Real>("initial_damage_random")
                                            : std::vector<Real>(_n_composition, 0.0)),
    _pf0(isParamValid("reference_fluid_pressure") ? getLynxParam<Real>("reference_fluid_pressure") : std::vector<Real>(_n_composition, 1.0)),
    _damage(declareADProperty<Real>("damage")),
    _damage_old(getMaterialPropertyOld<Real>("damage")),
    _damage_drive(declareADProperty<Real>("damage_force")),
    _plastic_strain_incr(declareADProperty<RankTwoTensor>("plastic_strain_increment")),
    _damage_incr(declareADProperty<Real>("damage_increment")),
    _yield_function(declareADProperty<Real>("plastic_yield_function")),
    _damage_poro_mech(declareADProperty<Real>("damage_poro_mechanical"))
{
}

void
LynxADDamageModelBase::setQp(unsigned int qp)
{
  _qp = qp;
}

void
LynxADDamageModelBase::initQpStatefulProperties()
{
  _damage[_qp] = averageProperty(_damage0) + (MooseRandom::rand() - 0.5) * averageProperty(_damage0_rand); 
}