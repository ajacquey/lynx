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

class LynxTemperatureCelsiusAux;

template <>
InputParameters validParams<LynxTemperatureCelsiusAux>();

class LynxTemperatureCelsiusAux : public AuxKernel
{
public:
  LynxTemperatureCelsiusAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _T_K;
};