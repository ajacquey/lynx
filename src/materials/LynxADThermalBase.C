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

#include "LynxADThermalBase.h"

InputParameters
LynxADThermalBase::validParams()
{
  InputParameters params = LynxADMaterialBase::validParams();
  params.addClassDescription("Base class for calculating the thermal properties.");
  params.addCoupledVar("porosity", 0.0, "The porosity auxiliary variable.");
  params.addParam<std::vector<Real>>("radiogenic_heat_production",
                                     "The radiogenic heat production value.");
  return params;
}

LynxADThermalBase::LynxADThermalBase(const InputParameters & parameters)
  : LynxADMaterialBase(parameters),
    // _coupled_porosity(isCoupled("porosity")),
    _porosity(coupledValue("porosity")),
    _heat_source(isParamValid("radiogenic_heat_production")
                     ? getLynxParam<Real>("radiogenic_heat_production")
                     : std::vector<Real>(_n_composition, 0.0)),
    _coupled_dens(hasMaterialProperty<RealVectorValue>("gravity_vector")),
    _rho_f(_coupled_dens ? &getADMaterialProperty<Real>("fluid_density") : nullptr),
    _rho_s(_coupled_dens ? &getADMaterialProperty<Real>("solid_density") : nullptr),
    _thermal_diff(declareADProperty<Real>("thermal_diffusivity")),
    _rhoC_b(declareADProperty<Real>("bulk_specific_heat")),
    _rhoC_f(declareADProperty<Real>("fluid_specific_heat")),
    _thermal_exp(declareProperty<Real>("thermal_expansion_coefficient")),
    _radiogenic_heat(declareProperty<Real>("radiogenic_heat_production")),
    _c_f(_fe_problem.getMaxQps()),
    _c_s(_fe_problem.getMaxQps()),
    _lambda_f(_fe_problem.getMaxQps()),
    _lambda_s(_fe_problem.getMaxQps()),
    _beta_f(_fe_problem.getMaxQps()),
    _beta_s(_fe_problem.getMaxQps())
{
}

void
LynxADThermalBase::computeQpProperties()
{
  computeQpSpecificHeat();
  computeQpThermalDiff();
  computeQpThermalExpansion();
  computeQpThermalSource();
}

void
LynxADThermalBase::computeQpSpecificHeat()
{
  computeQpHeatCap();

  _rhoC_f[_qp] = _coupled_dens ? (*_rho_f)[_qp] * _c_f[_qp] : 0.0;
  _rhoC_f[_qp] = 0.0;
  ADReal rhoC_s = _coupled_dens ? (*_rho_s)[_qp] * _c_s[_qp] : 0.0;
  _rhoC_b[_qp] = computeADMixtureProperty(_rhoC_f[_qp], rhoC_s);
}

void
LynxADThermalBase::computeQpThermalDiff()
{
  computeQpThermalCond();

  _thermal_diff[_qp] = computeMixtureProperty(_lambda_f[_qp], _lambda_s[_qp]);
  if (_rhoC_b[_qp] != 0.0)
    _thermal_diff[_qp] /= _rhoC_b[_qp];
}

void
LynxADThermalBase::computeQpThermalExpansion()
{
  computeQpThermalExp();

  _thermal_exp[_qp] = computeMixtureProperty(_beta_f[_qp], _beta_s[_qp]);
}

void
LynxADThermalBase::computeQpThermalSource()
{
  // In this material, we compute only radiogenic heat production
  // In LynxElasticRheology, we compute shear heating and adiabatic heating
  _radiogenic_heat[_qp] = averageProperty(_heat_source);
}

Real
LynxADThermalBase::computeMixtureProperty(const Real fluid_prop, const Real solid_prop)
{
  return _porosity[_qp] * fluid_prop + (1.0 - _porosity[_qp]) * solid_prop;
}

ADReal
LynxADThermalBase::computeADMixtureProperty(const ADReal fluid_prop, const ADReal solid_prop)
{
  return _porosity[_qp] * fluid_prop + (1.0 - _porosity[_qp]) * solid_prop;
}