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

#include "LynxThermalBase.h"

template <>
InputParameters
validParams<LynxThermalBase>()
{
  InputParameters params = validParams<LynxMaterialBase>();
  params.addClassDescription("Base class for calculating the thermal properties.");
  params.addRequiredCoupledVar("temperature", "The temperature variable.");
  params.addCoupledVar("porosity", "The porosity auxiliary variable.");
  params.addParam<std::vector<Real>>("radiogenic_heat_production",
                                     "The radiogenic heat production value.");
  return params;
}

LynxThermalBase::LynxThermalBase(const InputParameters & parameters)
  : LynxMaterialBase(parameters),
    _temp_dot(_fe_problem.isTransient() ? coupledDot("temperature") : _zero),
    _porosity(isCoupled("porosity") ? coupledValue("porosity") : _zero),
    _heat_source(isParamValid("radiogenic_heat_production")
                     ? getLynxParam<Real>("radiogenic_heat_production")
                     : std::vector<Real>(_n_composition, 0.0)),
    _rho_f(getDefaultMaterialProperty<Real>("fluid_density")),
    _rho_s(getDefaultMaterialProperty<Real>("solid_density")),
    _thermal_diff(declareProperty<Real>("thermal_diffusivity")),
    _rhoC_b(declareProperty<Real>("bulk_specific_heat")),
    _rhoC_f(declareProperty<Real>("fluid_specific_heat")),
    _thermal_exp(declareProperty<Real>("thermal_expansion_coefficient")),
    _thermal_strain_incr(declareProperty<RankTwoTensor>("thermal_strain_increment")),
    _dthermal_strain_dtemp(declareProperty<RankTwoTensor>("dthermal_strain_dtemp")),
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
LynxThermalBase::computeQpProperties()
{
  computeQpSpecificHeat();
  computeQpThermalDiff();
  computeQpThermalExpansion();
  computeQpThermalStrain();
  computeQpThermalSource();
}

void
LynxThermalBase::computeQpSpecificHeat()
{
  computeQpHeatCap();

  _rhoC_f[_qp] = _rho_f[_qp] * _c_f[_qp];
  _rhoC_f[_qp] = 0.0;
  _rhoC_b[_qp] = computeMixtureProperty(_rhoC_f[_qp], _rho_s[_qp] * _c_s[_qp]);
}

void
LynxThermalBase::computeQpThermalDiff()
{
  computeQpThermalCond();

  _thermal_diff[_qp] = computeMixtureProperty(_lambda_f[_qp], _lambda_s[_qp]);
  if (_rhoC_b[_qp] != 0.0)
    _thermal_diff[_qp] /= _rhoC_b[_qp];
}

void
LynxThermalBase::computeQpThermalExpansion()
{
  computeQpThermalExp();

  _thermal_exp[_qp] = computeMixtureProperty(_beta_f[_qp], _beta_s[_qp]);
}

void
LynxThermalBase::computeQpThermalStrain()
{
  _thermal_strain_incr[_qp] = RankTwoTensor();
  _thermal_strain_incr[_qp].addIa(_thermal_exp[_qp] * _temp_dot[_qp] * _dt / 3.0);

  _dthermal_strain_dtemp[_qp].zero();
  _dthermal_strain_dtemp[_qp].addIa(_thermal_exp[_qp] / 3.0);
  // + dthermal_exp_dtemp * _temp_dot / 3
}

void
LynxThermalBase::computeQpThermalSource()
{
  // In this material, we compute only radiogenic heat production
  // In LynxElasticRheology, we compute shear heating and adiabatic heating
  _radiogenic_heat[_qp] = averageProperty(_heat_source);
}

Real
LynxThermalBase::computeMixtureProperty(const Real fluid_prop, const Real solid_prop)
{
  return _porosity[_qp] * fluid_prop + (1.0 - _porosity[_qp]) * solid_prop;
}
