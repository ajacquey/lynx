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

#include "DirichletBCBase.h"

class LynxVelocityBC : public DirichletBCBase
{
public:
  static InputParameters validParams();
  LynxVelocityBC(const InputParameters & parameters);

protected:
  virtual Real computeQpValue() override;

  const VariableValue & _u_old;
  const Real & _value;
  const Function * _function;
};