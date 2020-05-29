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

#include "ADKernel.h"

class LynxADHeatSources : public ADKernel
{
public:
  static InputParameters validParams();
  LynxADHeatSources(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _coeff_Hs;
  unsigned int _nvel;
  std::vector<const VariableValue *> _vel;
  const ADVariableGradient & _grad_pressure;
  const ADMaterialProperty<Real> & _rhoC_b;
  const MaterialProperty<Real> & _radiogenic_heat;
  const bool _has_inelastic_heat_mat;
  const ADMaterialProperty<Real> * _inelastic_heat_mat;
  const bool _coupled_inelastic_heat;
  const VariableValue & _inelastic_heat;
  const MaterialProperty<Real> & _thermal_exp;
};