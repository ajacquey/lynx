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

class LynxDamageDeformation : public LynxDeformation
{
public:
  LynxDamageDeformation(const InputParameters & parameters);
  virtual ~LynxDamageDeformation() {}

protected:
  virtual void plasticCorrection(Real & pressure, RankTwoTensor & stress_dev);
  virtual void convexPlasticCorrection(Real & pressure, RankTwoTensor & stress_dev);
  virtual Real convexReferencePressure(const Real & p_tr, const Real & q_tr);
  virtual Real convexPlasticYield2(const Real & rho);
  virtual Real dConvexPlasticYield2(const Real & rho);
  virtual Real convexPlasticYield2(const Real & pressure, const Real & eqv_stress);
  virtual Real dConvexPlasticYield2_dp(const Real & pressure, const Real & eqv_stress);
  virtual Real dConvexPlasticYield2_dq(const Real & pressure, const Real & eqv_stress);
  virtual Real getConvexProjection(const Real & x1, const Real & x2);

  // Coupled variables
  bool _coupled_dam;
  const VariableValue & _damage;
  const VariableValue & _damage_old;
  bool _coupled_phi;
  const VariableValue & _porosity;

  // Strain parameters
  
  // Elastic moduli parameters
  const std::vector<Real> _damage_modulus;
  
  // Creep parameters

  // Plastic parameters
  const std::vector<Real> _critical_pressure;

  // Convex yield parameters
  Real _p_cr;
  Real _xi0;
  Real _p_r;
  Real _q_r;
  Real _dp_r_dp_tr;
  Real _dp_r_dq_tr;
  Real _p_tr;
  Real _q_tr;
  Real _rho_tr;
  Real _rp;
  Real _rq;
  Real _p_k;

  // Strain properties

  // Viscous properties

  // Plastic properties

  // Stress properties
};

#endif // LYNXDAMAGEDEFORMATION_H
