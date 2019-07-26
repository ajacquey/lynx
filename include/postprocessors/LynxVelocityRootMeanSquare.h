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

#include "ElementIntegralPostprocessor.h"

class LynxVelocityRootMeanSquare;

template <>
InputParameters validParams<LynxVelocityRootMeanSquare>();

class LynxVelocityRootMeanSquare : public ElementIntegralPostprocessor
{
public:
  LynxVelocityRootMeanSquare(const InputParameters & parameters);
  virtual Real getValue() override;

protected:
  virtual Real computeQpIntegral() override;

  unsigned int _n_vel;
  std::vector<const VariableValue *> _vel;
};