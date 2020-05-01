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

#include "ElementAverageValue.h"

class LynxElementAverageValue : public ElementAverageValue
{
public:
  static InputParameters validParams();
  LynxElementAverageValue(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  const VariableValue & _value_old;
  const VariableValue & _value_older;
};