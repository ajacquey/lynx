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

#include "AuxKernel.h"
#include "DerivativeMaterialInterface.h"

class LynxDarcyVelocityAux;

template <>
InputParameters validParams<LynxDarcyVelocityAux>();

class LynxDarcyVelocityAux : public DerivativeMaterialInterface<AuxKernel>
{
public:
  LynxDarcyVelocityAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const unsigned int _component;
  const VariableGradient & _grad_pf;
  const MaterialProperty<Real> & _fluid_mobility;
  const MaterialProperty<Real> & _rho_f;
  const MaterialProperty<RealVectorValue> & _gravity;
};