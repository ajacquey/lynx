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

#include "LynxADHeatSources.h"

registerMooseObject("LynxApp", LynxADHeatSources);

InputParameters
LynxADHeatSources::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription(
      "Heat generation by radiogenic, shear heating and adiabatic processes.");
  params.addParam<Real>(
      "coeff_shear_heating", 0.0, "The coefficient in front of the shear heating generation.");
  params.addCoupledVar(
      "inelastic_heat",
      "The auxiliary variable holding the inelastic heat value for running in a subApp.");
  params.addCoupledVar("velocities", "The string of velocities suitable for the problem statement");
  params.addCoupledVar("pressure", "The pressure variable");
  return params;
}

LynxADHeatSources::LynxADHeatSources(const InputParameters & parameters)
  : ADKernel(parameters),
    _coeff_Hs(getParam<Real>("coeff_shear_heating")),
    _nvel(coupledComponents("velocities")),
    _vel(3),
    _grad_pressure(isCoupled("pressure") ? adCoupledGradient("pressure") : _ad_grad_zero),
    _rhoC_b(getADMaterialProperty<Real>("bulk_specific_heat")),
    _radiogenic_heat(getMaterialProperty<Real>("radiogenic_heat_production")),
    _has_inelastic_heat_mat(hasADMaterialProperty<Real>("inelastic_heat")),
    _inelastic_heat_mat(_has_inelastic_heat_mat ? &getADMaterialProperty<Real>("inelastic_heat")
                                                : nullptr),
    _coupled_inelastic_heat(isCoupled("inelastic_heat")),
    _inelastic_heat(_coupled_inelastic_heat ? coupledValue("inelastic_heat") : _zero),
    _thermal_exp(getMaterialProperty<Real>("thermal_expansion_coefficient"))
{
  for (unsigned i = 0; i < _nvel; ++i)
    _vel[i] = &coupledValue("velocities", i);
  for (unsigned i = _nvel; i < 3; ++i)
    _vel[i] = &_zero;
}

ADReal
LynxADHeatSources::computeQpResidual()
{
  Real Hr = _radiogenic_heat[_qp];
  ADReal Hs = _coeff_Hs;
  if (_coupled_inelastic_heat)
    Hs *= _inelastic_heat[_qp];
  else if (_has_inelastic_heat_mat)
    Hs *= (*_inelastic_heat_mat)[_qp];
  else
    Hs *= 0.0;
  RealVectorValue vel((*_vel[0])[_qp], (*_vel[1])[_qp], (*_vel[2])[_qp]);
  ADReal Ha = _thermal_exp[_qp] * _u[_qp] * vel * _grad_pressure[_qp];

  ADReal heat_sources = Hr + Hs + Ha;
  if (_rhoC_b[_qp] != 0.0)
    heat_sources /= _rhoC_b[_qp];

  return -heat_sources * _test[_i][_qp];
}