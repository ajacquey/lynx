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

#include "LynxDensityConstant.h"

registerMooseObject("LynxApp", LynxDensityConstant);

template <>
InputParameters
validParams<LynxDensityConstant>()
{
  InputParameters params = validParams<LynxDensityBase>();
  params.addClassDescription("Material calculating densities as constant values.");
  return params;
}

LynxDensityConstant::LynxDensityConstant(const InputParameters & parameters)
  : LynxDensityBase(parameters)
{
}

void
LynxDensityConstant::computeQpProperties()
{
  computeQpGravity();
  _rho_f[_qp] = averageProperty(_fluid_density);
  _rho_s[_qp] = averageProperty(_solid_density);
  _rho_b[_qp] = _porosity[_qp] * _rho_f[_qp] + (1.0 - _porosity[_qp]) * _rho_s[_qp];
  _reference_rho_b[_qp] = _rho_b[_qp];
}
