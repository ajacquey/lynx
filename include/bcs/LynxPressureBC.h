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

#ifndef LYNXPRESSUREBC_H
#define LYNXPRESSUREBC_H

#include "IntegratedBC.h"

class LynxPressureBC;
class Function;

template <>
InputParameters validParams<LynxPressureBC>();

class LynxPressureBC : public IntegratedBC
{
public:
  LynxPressureBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  const int _component;
  const Real _value;
  const Function * _function;
};

#endif // LYNXPRESSUREBC_H
