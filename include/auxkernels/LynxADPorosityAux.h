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

class LynxADPorosityAux : public AuxKernel
{
public:
  static InputParameters validParams();
  LynxADPorosityAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  virtual Real computeEvDot();
  virtual Real computeEvInDot();

  bool _coupled_pf;
  const VariableValue & _pf_dot;
  const MaterialProperty<Real> & _biot;
  const MaterialProperty<Real> & _C_d;
  const ADMaterialProperty<RankTwoTensor> & _strain_increment;
  const bool _has_viscous;
  const ADMaterialProperty<RankTwoTensor> * _viscous_strain_incr;
  const bool _has_plastic;
  const ADMaterialProperty<RankTwoTensor> * _plastic_strain_incr;
  const bool _has_damage;
  const ADMaterialProperty<Real> * _damage_poro_mech;
};