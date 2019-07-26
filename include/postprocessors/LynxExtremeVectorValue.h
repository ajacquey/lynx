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
#ifndef LYNXEXTREMEVECTORVALUE_H
#define LYNXEXTREMEVECTORVALUE_H

#include "ElementExtremeValue.h"
#include "DerivativeMaterialInterface.h"

class LynxExtremeVectorValue;

template <>
InputParameters validParams<LynxExtremeVectorValue>();

class LynxExtremeVectorValue : public DerivativeMaterialInterface<ElementExtremeValue>
{
public:
  LynxExtremeVectorValue(const InputParameters & parameters);

protected:
  virtual void computeQpValue() override;
  const VariableValue & _v;
  const VariableValue & _w;
};

#endif // LYNXEXTREMESCALARVALUE_H
