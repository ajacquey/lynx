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

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

class LynxHydroDarcy;

template <>
InputParameters validParams<LynxHydroDarcy>();

class LynxHydroDarcy : public DerivativeMaterialInterface<Kernel>
{
public:
  LynxHydroDarcy(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  // virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const MaterialProperty<Real> & _fluid_mobility;
  const MaterialProperty<RealVectorValue> & _gravity;
  const MaterialProperty<Real> & _rho_f;
};