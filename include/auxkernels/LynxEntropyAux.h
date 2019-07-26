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

class LynxEntropyAux;
template <>
InputParameters validParams<LynxEntropyAux>();

class LynxEntropyAux : public DerivativeMaterialInterface<AuxKernel>
{
public:
  LynxEntropyAux(const InputParameters & parameters);
  virtual ~LynxEntropyAux() {}

protected:
  virtual Real computeValue();
  const VariableValue & _var_old;
  const VariableValue & _var_older;

  const PostprocessorValue & _pp_max_var;
  const PostprocessorValue & _pp_min_var;
};