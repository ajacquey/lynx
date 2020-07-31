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

#include "LynxADMaterialBase.h"

class LynxADDamageModelBase : public LynxADMaterialBase
{
public:
  static InputParameters validParams();
  LynxADDamageModelBase(const InputParameters & parameters);
  void setQp(unsigned int qp);
  // virtual void elasticGuess(ADRankTwoTensor & stress,
  //                           const ADRankTwoTensor & stress_old,
  //                           const ADRankFourTensor & Cijkl,
  //                           const ADRankTwoTensor & elastic_strain_old,
  //                           const ADRankTwoTensor & elastic_strain_incr) = 0;
  // virtual void damageUpdate(ADRankTwoTensor & stress, ADRankTwoTensor & elastic_strain_incr) = 0;
  virtual void damageUpdate(ADRankTwoTensor & stress,
                            const ADRankFourTensor & Cijkl,
                            const ADRankTwoTensor & elastic_strain_old,
                            ADRankTwoTensor & elastic_strain_incr) = 0;
  virtual void updateElasticModuli(ADReal & K, ADReal & G, const ADRankTwoTensor & elastic_strain) = 0;
  void resetQpProperties() final {}
  void resetProperties() final {}

protected:
  virtual void initQpStatefulProperties() override;

  const ADVariableValue & _pf;
  const Real _abs_tol;
  const Real _rel_tol;
  const unsigned int _max_its;
  const std::vector<Real> _damage0;
  const std::vector<Real> _damage0_rand;
  const std::vector<Real> _pf0;

  // Damage properties
  ADMaterialProperty<Real> & _damage;
  const MaterialProperty<Real> & _damage_old;
  ADMaterialProperty<Real> & _damage_drive;
  ADMaterialProperty<RankTwoTensor> & _plastic_strain_incr;
  ADMaterialProperty<Real> & _damage_incr;
  ADMaterialProperty<Real> & _yield_function;
  ADMaterialProperty<Real> & _damage_poro_mech;
};