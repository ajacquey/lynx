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

#ifndef LYNXEXTRAPOLATEVALUE_H
#define LYNXEXTRAPOLATEVALUE_H

#include "ElementExtremeValue.h"
#include "DerivativeMaterialInterface.h"

class LynxExtrapolateValue;

template <>
InputParameters validParams<LynxExtrapolateValue>();

class LynxExtrapolateValue : public DerivativeMaterialInterface<ElementExtremeValue>
{
public:
  LynxExtrapolateValue(const InputParameters & parameters);

protected:
  virtual void computeQpValue() override;

  const VariableValue & _value_old;
  const VariableValue & _value_older;
};

#endif // LYNXEXTRAPOLATEVALUE_H
