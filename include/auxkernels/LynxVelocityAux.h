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

#ifndef LYNXVELOCITYAUX_H
#define LYNXVELOCITYAUX_H

#include "AuxKernel.h"

class LynxVelocityAux;

template <>
InputParameters validParams<LynxVelocityAux>();

class LynxVelocityAux : public AuxKernel
{
public:
  LynxVelocityAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _disp;
  const VariableValue & _disp_old;
};

#endif // LYNXVELOCITYAUX_H
