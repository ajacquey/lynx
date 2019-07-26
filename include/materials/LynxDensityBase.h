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
