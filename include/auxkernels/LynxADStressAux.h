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

#include "LynxADStressAuxBase.h"

class LynxADStressAux : public LynxADStressAuxBase
{
public:
  static InputParameters validParams();
  LynxADStressAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:
  const bool _coupled_pt;
  const ADVariableValue & _pt;
  const unsigned int _i;
  const unsigned int _j;
};