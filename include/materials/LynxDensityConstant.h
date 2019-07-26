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

#include "LynxDensityBase.h"

class LynxDensityConstant;

template <>
InputParameters validParams<LynxDensityConstant>();

class LynxDensityConstant : public LynxDensityBase
{
public:
  LynxDensityConstant(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
};