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

#include "LynxDeformation.h"

class LynxDamageDeformation;

template <>
InputParameters validParams<LynxDamageDeformation>();

class LynxDamageDeformation : public LynxDeformationBase
{
public:
  LynxDamageDeformation(const InputParameters & parameters);
  virtual ~LynxDamageDeformation() {}

protected:
  virtual void initializeQpDeformation() override;
  virtual void plasticCorrection(Real & pressure, RankTwoTensor & stress_dev) override;
  virtual void damageCorrection() override;
  virtual RankFourTensor damageTangentOperator(const RankFourTensor & tme) override;
  virtual Real computePlasticityYield(const Real & pressure, const Real & eqv_stress);
  virtual Real plasticIncrement(const Real & /*pressure*/, const Real & eqv_stress);
  virtual Real computeConvexPlasticityYield2(const Real & pressure, const Real & eqv_stress);
  virtual Real computeConvexPlasticityYield2(const Real & rho);
  virtual Real convexPlasticIncrement(Real & vol_plastic_incr, Real & eqv_plastic_incr);
  virtual void computeDamageProperties(const Real & pressure, const Real & eqv_stress);
  virtual void updateDamageParameters();
  virtual void updateDamageConvexParameters(const Real & pressure, const Real & eqv_stress);
  virtual Real convexReferencePressure();
  virtual Real dConvexPlasticYield2(const Real & rho);
  virtual Real dConvexPlasticYield2_dp(const Real & pressure, const Real & eqv_stress);
  virtual Real dConvexPlasticYield2_dq(const Real & pressure, const Real & eqv_stress);
  virtual Real getConvexProjection(const Real & x1, const Real & x2);
  virtual Real strainRatio(const RankTwoTensor & elastic_strain);
  virtual RankTwoTensor rotatedElasticStrain(const RankTwoTensor & elastic_strain);
  virtual void computeQpThermalSources() override;

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
  Real _dyield_dp_tr;
  Real _dyield_dq_tr;

  // Stress properties
  MaterialProperty<RankTwoTensor> & _dstress_ddamage;
  MaterialProperty<RankTwoTensor> & _ddamage_rate_dstrain;

  // Damage properties
  MaterialProperty<Real> & _damage_rate;
  MaterialProperty<Real> & _damage_heat;
};