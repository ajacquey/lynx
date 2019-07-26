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

#include "IntegratedBC.h"

class LynxHeatFluxBC;
class Function;
class MooseRandom;

template <>
InputParameters validParams<LynxHeatFluxBC>();

class LynxHeatFluxBC : public IntegratedBC
{
public:
  LynxHeatFluxBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  const Real & _value;
  const Function * _function;
  Real _rand_per;
  const MaterialProperty<Real> & _rhoC_b;
};