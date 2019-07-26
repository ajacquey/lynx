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

#include "LynxHeatFluxAux.h"

registerMooseObject("LynxApp", LynxHeatFluxAux);

template <>
InputParameters
validParams<LynxHeatFluxAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Calculates the heat flux in each element for the given direction.");
  params.addRequiredCoupledVar("temperature", "The temperature variable.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "component",
      "component >= 0 & component <= 2",
      "Integer corresponding to direction of the heat flux (0, 1, 2.");
  return params;
}

LynxHeatFluxAux::LynxHeatFluxAux(const InputParameters & parameters)
  : DerivativeMaterialInterface<AuxKernel>(parameters),
    _component(getParam<unsigned int>("component")),
    _grad_T(coupledGradient("temperature")),
    _rhoC_b(getDefaultMaterialProperty<Real>("bulk_specific_heat")),
    _thermal_diff(getDefaultMaterialProperty<Real>("thermal_diffusivity"))
{
}

Real
LynxHeatFluxAux::computeValue()
{
  if (_rhoC_b[_qp] == 0.0)
    _cond_bulk = _thermal_diff[_qp];
  else
    _cond_bulk = _thermal_diff[_qp] * _rhoC_b[_qp];

  return -1 * _cond_bulk * _grad_T[_qp](_component);
}
