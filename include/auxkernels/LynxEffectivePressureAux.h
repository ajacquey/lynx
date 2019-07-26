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

#ifndef LYNXEFFECTIVEPRESSUREAUX_H
#define LYNXEFFECTIVEPRESSUREAUX_H

#include "LynxStressAuxBase.h"

class LynxEffectivePressureAux;

template <>
InputParameters validParams<LynxEffectivePressureAux>();

class LynxEffectivePressureAux : public LynxStressAuxBase
{
public:
  LynxEffectivePressureAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
};

#endif // LYNXEFFECTIVEPRESSUREAUX_H
