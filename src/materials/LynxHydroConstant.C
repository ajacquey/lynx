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

#include "LynxHydroConstant.h"

registerMooseObject("LynxApp", LynxHydroConstant);

template <>
InputParameters
validParams<LynxHydroConstant>()
{
  InputParameters params = validParams<LynxHydroBase>();
  params.addClassDescription("Constant thermal properties.");
  params.addRequiredParam<std::vector<Real>>("permeability", "The permeability of the matrix.");
  params.addRequiredParam<std::vector<Real>>("fluid_viscosity", "The viscosity of the fluid.");
  params.addParam<std::vector<Real>>("fluid_modulus", "The bulk modulus of the fluid phase.");
  params.addParam<std::vector<Real>>("solid_modulus", "The bulk modulus of the solid phase.");
  return params;
}

LynxHydroConstant::LynxHydroConstant(const InputParameters & parameters)
  : LynxHydroBase(parameters),
    _perm(getLynxParam<Real>("permeability")),
    _fluid_viscosity(getLynxParam<Real>("fluid_viscosity"))
{
  if (isParamValid("fluid_modulus"))
  {
    _fluid_compr = std::vector<Real>(_n_composition, 0.0);
    std::vector<Real> fluid_modulus = getLynxParam<Real>("fluid_modulus");
    for (unsigned int i = 0; i < _n_composition; ++i)
      if (fluid_modulus[i] != 0.0)
        _fluid_compr[i] = 1.0 / fluid_modulus[i];
  }
  else
    _fluid_compr = std::vector<Real>(_n_composition, 0.0);

  if (isParamValid("solid_modulus"))
  {
    _solid_compr = std::vector<Real>(_n_composition, 0.0);
    std::vector<Real> solid_modulus = getLynxParam<Real>("solid_modulus");
    for (unsigned int i = 0; i < _n_composition; ++i)
      if (solid_modulus[i] != 0.0)
        _solid_compr[i] = 1.0 / solid_modulus[i];
  }
  else
    _solid_compr = std::vector<Real>(_n_composition, 0.0);
}

void
LynxHydroConstant::computeQpFluidCompressibility()
{
  _C_f[_qp] = averageProperty(_fluid_compr);
}

void
LynxHydroConstant::computeQpSolidCompressibility()
{
  _C_s[_qp] = averageProperty(_solid_compr);
}

void
LynxHydroConstant::computeQpPermeability()
{
  _k[_qp] = averageProperty(_perm);
}

void
LynxHydroConstant::computeQpFluidViscosity()
{
  _eta_f[_qp] = averageProperty(_fluid_viscosity);
}