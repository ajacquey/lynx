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

#pragma once

#include "LynxADMaterialBase.h"

class LynxADDensityBase : public LynxADMaterialBase
{
public:
  static InputParameters validParams();
  LynxADDensityBase(const InputParameters & parameters);

protected:
  virtual void computeQpGravity();

  const ADVariableValue & _porosity;

  bool _has_gravity;
  Real _g;
  const std::vector<Real> _fluid_density;
  const std::vector<Real> _solid_density;

  ADMaterialProperty<RealVectorValue> & _gravity;
  ADMaterialProperty<Real> & _rho_f;
  ADMaterialProperty<Real> & _rho_s;
  ADMaterialProperty<Real> & _rho_b;
  ADMaterialProperty<Real> & _reference_rho_b;
};