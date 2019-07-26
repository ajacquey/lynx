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

#ifndef LYNXCOMBOPHASESAUX_H
#define LYNXCOMBOPHASESAUX_H

#include "AuxKernel.h"

class LynxComboPhasesAux;
template <>
InputParameters validParams<LynxComboPhasesAux>();

class LynxComboPhasesAux : public AuxKernel
{
public:
  LynxComboPhasesAux(const InputParameters & parameters);
  virtual ~LynxComboPhasesAux() {}

protected:
  virtual Real computeValue();

  unsigned _n_composition;
  std::vector<const VariableValue *> _compositional_phases;
};

#endif // LYNXCOMBOPHASESAUX_H
