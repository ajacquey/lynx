/******************************************************************************/
/*                       LYNX, a MOOSE-based application                      */
/*                                                                            */
/*          Copyright (C) 2017 by Antoine B. Jacquey and Mauro Cacace         */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*    This program is free software: you can redistribute it and/or modify    */
/*    it under the terms of the GNU General Public License as published by    */
/*      the Free Software Foundation, either version 3 of the License, or     */
/*                     (at your option) any later version.                    */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program. If not, see <http://www.gnu.org/licenses/>     */
/******************************************************************************/

#ifndef LYNXDAMAGEDEFORMATIONNEW_H
#define LYNXDAMAGEDEFORMATIONNEW_H

#include "LynxDeformation.h"

class LynxDamageDeformationNew;

template <>
InputParameters validParams<LynxDamageDeformationNew>();

class LynxDamageDeformationNew : public LynxDeformationBase
{
public:
  LynxDamageDeformationNew(const InputParameters & parameters);
  virtual ~LynxDamageDeformationNew() {}

protected:
  virtual void initQpStatefulProperties() override;
  virtual void initializeQpDeformation() override;
  // virtual void elasticModuli() override;
  virtual void plasticCorrection(Real & pressure, RankTwoTensor & stress_dev) override;
  virtual void damageCorrection() override;
  virtual RankFourTensor damageTangentOperator(const RankTwoTensor & flow_direction,
                                               const RankFourTensor & tme) override;
  virtual Real computePlasticityYield(const Real & pressure, const Real & eqv_stress);
  virtual Real plasticIncrement(const Real & /*pressure*/, const Real & eqv_stress);
  virtual Real computeConvexPlasticityYield2(const Real & pressure, const Real & eqv_stress);
  virtual Real computeConvexPlasticityYield2(const Real & rho);
  virtual Real convexPlasticIncrement(Real & vol_plastic_incr, Real & eqv_plastic_incr);
  virtual void computeDamageProperties(const Real & pressure, const Real & eqv_stress);
  virtual void updateDamageParameters();
  virtual void initializeDamageParameters();
  virtual void updateDamageConvexParameters(const Real & pressure, const Real & eqv_stress);
  virtual Real convexReferencePressure();
  virtual Real dConvexPlasticYield2(const Real & rho);
  virtual Real dConvexPlasticYield2_dp(const Real & pressure, const Real & eqv_stress);
  virtual Real dConvexPlasticYield2_dq(const Real & pressure, const Real & eqv_stress);
  virtual Real getConvexProjection(const Real & x1, const Real & x2);
  virtual Real strainRatio(const RankTwoTensor & elastic_strain);
  virtual RankTwoTensor rotatedElasticStrain(const RankTwoTensor & elastic_strain);

  // Coupled variables
  bool _coupled_dam;
  const VariableValue & _damage;
  const VariableValue & _damage_old;
  bool _coupled_phi;
  const VariableValue & _porosity;

  // Elastic moduli parameters
  const std::vector<Real> _damage_modulus;

  // Plastic parameters
  const std::vector<Real> _friction_angle;
  const std::vector<Real> _cohesion;
  const std::vector<Real> _porous_coeff;
  const std::vector<Real> _porous_coeff_linear;
  std::vector<Real> _one_on_plastic_eta;
  std::vector<Real> _one_on_damage_eta;

  // Damage-Plasticity structure
  damage_plasticity * _damage_plasticity;

  // Damage-Plasticity utilities
  // Real _dyield_dp_tr;
  // Real _dyield_dq_tr;
  // RankTwoTensor _damaged_stress;
  // RankFourTensor _damaged_tensor;

  // Strain properties
  MaterialProperty<RankTwoTensor> & _elastic_strain;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;

  // Stress properties
  // MaterialProperty<RankTwoTensor> & _dstress_ddamage;

  // Damage properties
  MaterialProperty<Real> & _damage_rate;
};

#endif // LYNXDAMAGEDEFORMATIONNEW_H
