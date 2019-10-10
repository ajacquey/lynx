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

template <ComputeStage>
class LynxADPlasticModel;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;

declareADValidParams(LynxADPlasticModel);

template <ComputeStage compute_stage>
class LynxADPlasticModel : public LynxADMaterialBase<compute_stage>
{
public:
  LynxADPlasticModel(const InputParameters & parameters);
  void setQp(unsigned int qp);
  virtual void plasticUpdate(ADRankTwoTensor & stress_dev,
                             ADReal & pressure,
                             const ADReal & G,
                             const ADReal & K,
                             ADRankTwoTensor & elastic_strain_incr);
  void resetQpProperties() final {}
  void resetProperties() final {}

protected:
  virtual void initQpStatefulProperties() override;
  virtual ADReal plasticIncrement(const ADReal & eqv_stress, const ADReal & pressure, const ADReal G, const ADReal K);
  virtual void initPlasticParameters(const ADReal & pressure, const ADReal & K);
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
  ADMaterialProperty(Real) & _plastic_yield_function;
  ADMaterialProperty(RankTwoTensor) & _plastic_strain_incr;
  ADMaterialProperty(Real) * _intnl;
  const MaterialProperty<Real> * _intnl_old;

  // Plastic effective parameters
  ADReal _alpha_0;
  ADReal _alpha_res;
  ADReal _k_0;
  ADReal _k_res;
  ADReal _beta;
  ADReal _intnl_0;
  ADReal _intnl_lim;
  ADReal _one_on_eta;
  ADReal _alpha;
  ADReal _k;
  ADReal _H;

  usingMaterialBaseMembers;
};