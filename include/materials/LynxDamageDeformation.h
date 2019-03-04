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

#ifndef LYNXDAMAGEDEFORMATION_H
#define LYNXDAMAGEDEFORMATION_H

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
  virtual void initQpStatefulProperties() override;
  virtual void elasticModuli() override;
  virtual void plasticCorrection(Real & pressure, RankTwoTensor & stress_dev) override;
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
  const std::vector<Real> _porous_coeff;
  const std::vector<Real> _porous_coeff_linear;
  std::vector<Real> _one_on_plastic_eta;
  std::vector<Real> _one_on_damage_eta;

  // Damage-Plasticity structure
  damage_plasticity * _damage_plasticity;

  // Strain properties
  MaterialProperty<RankTwoTensor> & _elastic_strain;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
};

#endif // LYNXDAMAGEDEFORMATION_H
