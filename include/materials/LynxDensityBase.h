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

#ifndef LYNXDENSITYBASE_H
#define LYNXDENSITYBASE_H

#include "LynxMaterialBase.h"

class LynxDensityBase;

template <>
InputParameters validParams<LynxDensityBase>();

class LynxDensityBase : public LynxMaterialBase
{
public:
  LynxDensityBase(const InputParameters & parameters);
  virtual ~LynxDensityBase() {}

protected:
  virtual void computeQpGravity();

  const VariableValue & _porosity;

  bool _has_gravity;
  Real _g;
  const std::vector<Real> _fluid_density;
  const std::vector<Real> _solid_density;

  MaterialProperty<RealVectorValue> & _gravity;
  MaterialProperty<Real> & _rho_f;
  MaterialProperty<Real> & _rho_s;
  MaterialProperty<Real> & _rho_b;
  MaterialProperty<Real> & _reference_rho_b;
};

#endif // LYNXDENSITYBASE_H
