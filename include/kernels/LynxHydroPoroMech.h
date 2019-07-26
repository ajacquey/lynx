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

class LynxHydroPoroMech;

template <>
InputParameters validParams<LynxHydroPoroMech>();

class LynxHydroPoroMech : public DerivativeMaterialInterface<Kernel>
{
public:
  LynxHydroPoroMech(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  unsigned int _ndisp;
  const MaterialProperty<Real> & _poro_mech;
  const MaterialProperty<Real> & _poro_mech_jac;
  std::vector<unsigned int> _disp_var;
};