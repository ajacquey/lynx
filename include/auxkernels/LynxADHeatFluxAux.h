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

#include "AuxKernel.h"

class LynxADHeatFluxAux : public AuxKernel
{
public:
  static InputParameters validParams();
  LynxADHeatFluxAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const unsigned int _component;
  const VariableGradient & _grad_T;
  const ADMaterialProperty<Real> & _rhoC_b;
  const ADMaterialProperty<Real> & _thermal_diff;
  Real _cond_bulk;
};