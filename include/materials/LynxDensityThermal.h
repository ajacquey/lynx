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
