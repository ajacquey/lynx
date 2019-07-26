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

class LynxHoldStressBC;
class Function;

template <>
InputParameters validParams<LynxHoldStressBC>();

class LynxHoldStressBC : public IntegratedBC
{
public:
  LynxHoldStressBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  bool _coupled_pf;
  const VariableValue & _pf;
  const int _component;
  const Real _bulk_modulus;
  const Real _shear_modulus;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  const MaterialProperty<Real> & _biot_coeff;
};