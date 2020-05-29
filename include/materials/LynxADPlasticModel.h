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

class LynxADPlasticModel : public LynxADMaterialBase
{
public:
  static InputParameters validParams();
  LynxADPlasticModel(const InputParameters & parameters);
  void setQp(unsigned int qp);
  virtual void plasticUpdate(ADRankTwoTensor & stress_dev,
                             ADReal & pressure,
                             const Real & G,
                             const Real & K,
                             ADRankTwoTensor & elastic_strain_incr);
  void resetQpProperties() final {}
  void resetProperties() final {}

protected:
  virtual void initQpStatefulProperties() override;
  virtual ADReal plasticIncrement(const ADReal & eqv_stress, const ADReal & pressure, const Real G, const Real K);
  virtual void initPlasticParameters(const ADReal & pressure, const Real & K);
  virtual ADReal plasticYieldFunction(const ADReal & eqv_stress, const ADReal & pressure);

  // Plastic parameters
  const std::vector<Real> _friction_angle_0;
  const std::vector<Real> _cohesion_0;
  const std::vector<Real> _friction_angle_res;
  const std::vector<Real> _cohesion_res;
  const std::vector<Real> _dilation_angle;
  const std::vector<Real> _intnl_param_0;
  const std::vector<Real> _intnl_param_lim;
  std::vector<Real> _one_on_plastic_eta;
  const bool _has_hardening;

  // Plastic Properties
  ADMaterialProperty<Real> & _plastic_yield_function;
  ADMaterialProperty<RankTwoTensor> & _plastic_strain_incr;
  ADMaterialProperty<Real> * _intnl;
  const MaterialProperty<Real> * _intnl_old;

  // Plastic effective parameters
  Real _alpha_0;
  Real _alpha_res;
  Real _k_0;
  Real _k_res;
  Real _beta;
  Real _intnl_0;
  Real _intnl_lim;
  Real _one_on_eta;
  ADReal _alpha;
  ADReal _k;
  ADReal _H;
};