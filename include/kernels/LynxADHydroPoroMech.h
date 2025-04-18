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

#include "ADKernel.h"

class LynxADHydroPoroMech : public ADKernel
{
public:
  static InputParameters validParams();
  LynxADHydroPoroMech(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty<Real> & _poro_mech;
  const ADMaterialProperty<Real> & _biot;
  const bool _has_damage;
  const ADMaterialProperty<Real> * _damage_poro_mech;
};