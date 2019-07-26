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
#include "DerivativeMaterialInterface.h"

class LynxVariableRateAux;

template <>
InputParameters validParams<LynxVariableRateAux>();

class LynxVariableRateAux : public DerivativeMaterialInterface<AuxKernel>
{
public:
  LynxVariableRateAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _cvar;
  const VariableValue & _cvar_old;
  Real _tscale;
  bool _relative;
};