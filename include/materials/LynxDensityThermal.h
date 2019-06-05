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

#ifndef LYNXDENSITYTHERMAL_H
#define LYNXDENSITYTHERMAL_H

#include "LynxDensityBase.h"
#include "Function.h"

class LynxDensityThermal;

template <>
InputParameters validParams<LynxDensityThermal>();

class LynxDensityThermal : public LynxDensityBase
{
public:
  LynxDensityThermal(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const VariableValue & _temperature;

  MaterialProperty<Real> & _drho_dtemp;
  MaterialProperty<Real> & _dinvrho_dtemp;
  const std::vector<Real> _beta_fluid;
  const std::vector<Real> _beta_solid;
  Real _temp_ref;
  const Function * _temp_ref_fct;
};

#endif // LYNXDENSITYTHERMAL_H
