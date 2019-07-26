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

#ifndef LYNXWINKLERBC_H
#define LYNXWINKLERBC_H

#include "IntegratedBC.h"
#include "DerivativeMaterialInterface.h"

class LynxWinklerBC;
class Function;

template <>
InputParameters validParams<LynxWinklerBC>();

class LynxWinklerBC : public DerivativeMaterialInterface<IntegratedBC>
{
public:
  LynxWinklerBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  unsigned int _ndisp;
  std::vector<const VariableValue *> _disp;
  const int _component;
  const Real _value;
  const Function * _function;
  const Real _rho_ext;
  const Real _g;
  const MaterialProperty<Real> & _rho_b;
};

#endif // LYNXWINKLERBC_H
