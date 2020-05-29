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

class LynxADDarcyVelocityAux : public AuxKernel
{
public:
  static InputParameters validParams();
  LynxADDarcyVelocityAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const unsigned int _component;
  const VariableGradient & _grad_pf;
  const MaterialProperty<Real> & _fluid_mobility;
  const bool _coupled_grav;
  const MaterialProperty<RealVectorValue> * _gravity;
  const ADMaterialProperty<Real> * _rho_f;
};