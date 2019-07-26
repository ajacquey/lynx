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
#include "LynxDensityCompressible.h"

registerMooseObject("LynxApp", LynxDensityCompressible);

template <>
InputParameters
validParams<LynxDensityCompressible>()
{
  InputParameters params = validParams<LynxDensityBase>();
  params.addClassDescription("Material calculating compressible form of the solid density.");
  params.addCoupledVar("temperature", "The actual temperature variable.");
  params.addCoupledVar("reference_temperature", "The reference temperature variable.");
  params.addCoupledVar("lithostatic_pressure", "The lithostatic pressure.");
  params.addParam<std::vector<Real>>("beta_solid", "The solid thermal expansion coefficient.");
  // params.addParam<bool>(
  //     "temperature_from_multiapp", false, "Is the temperature obtained from a multiapp
  //     staging?");
  return params;
}

LynxDensityCompressible::LynxDensityCompressible(const InputParameters & parameters)
  : LynxDensityBase(parameters),
    _coupled_temp(isCoupled("temperature")),
    _temp(_coupled_temp ? coupledValue("temperature") : _zero),
    _reference_temperature(_coupled_temp ? coupledValue("reference_temperature") : _zero),
    _beta_solid(isParamValid("beta_solid") ? getLynxParam<Real>("beta_solid")
                                           : std::vector<Real>(_n_composition, 0.0)),
    // _temperature_from_multiapp(getParam<bool>("temperature_from_multiapp")),
    _drho_dtemp(declareProperty<Real>("drho_dtemp")),
    _drho_dev(declareProperty<Real>("drho_dev")),
    _dinvrho_dtemp(declareProperty<Real>("dinvrho_dtemp")),
    _dinvrho_dev(declareProperty<Real>("dinvrho_dev")),
    _coupled_plith(isCoupled("lithostatic_pressure")),
    _plith(_coupled_plith ? coupledValue("lithostatic_pressure") : _zero),
    _stress(getDefaultMaterialProperty<RankTwoTensor>("stress")),
    _K(getDefaultMaterialProperty<Real>("bulk_modulus"))
{
}

void
LynxDensityCompressible::computeQpProperties()
{
  computeQpGravity();

  // here calculate the density for the lithostatic pressure kernel
  _reference_rho_b[_qp] = averageProperty(_solid_density);

  // here calculate the density for the dynamic internal component
  _rho_s[_qp] = averageProperty(_solid_density);

  // augment the density for the volumetric deformation contribution
  Real one_on_K = _K[_qp] != 0.0 ? 1.0 / _K[_qp] : 0.0;
  Real dynamic_pressure = -_stress[_qp].trace() / 3.0 - _plith[_qp];
  _rho_s[_qp] += averageProperty(_solid_density) * one_on_K * dynamic_pressure;

  // augment the density for the thermal contribution
  Real dT = _temp[_qp] - _reference_temperature[_qp];
  _rho_s[_qp] -= averageProperty(_solid_density) * averageProperty(_beta_solid) * dT;
  // update the bulk density
  _rho_b[_qp] = _rho_s[_qp];

  // if (!_temperature_from_multiapp)
  // {
  Real drho_s_dtemp = -averageProperty(_solid_density) * averageProperty(_beta_solid);
  Real drho_s_dev = averageProperty(_solid_density);
  _drho_dtemp[_qp] = drho_s_dtemp;
  _drho_dev[_qp] = drho_s_dev;
  _dinvrho_dtemp[_qp] =
      (_rho_b[_qp] != 0.0) ? -1.0 / Utility::pow<2>(_rho_b[_qp]) * _drho_dtemp[_qp] : 0.0;
  _dinvrho_dev[_qp] =
      (_rho_b[_qp] != 0.0) ? -1.0 / Utility::pow<2>(_rho_b[_qp]) * _drho_dev[_qp] : 0.0;
  // }
}
