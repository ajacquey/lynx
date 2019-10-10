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
class LynxADCreepModel;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;

declareADValidParams(LynxADCreepModel);

template <ComputeStage compute_stage>
class LynxADCreepModel : public LynxADMaterialBase<compute_stage>
{
public:
  LynxADCreepModel(const InputParameters & parameters);
  void setQp(unsigned int qp);
  virtual void creepUpdate(ADRankTwoTensor & stress_dev,
                           const ADReal & pressure,
                           const ADReal & G,
                           ADRankTwoTensor & elastic_strain_incr);
  void resetQpProperties() final {}
  void resetProperties() final {}

protected:
  virtual ADReal viscousIncrement(const ADReal & pressure, const ADReal & eqv_stress, const ADReal & G);
  virtual void initCreepParameters(const ADReal & pressure);
  virtual ADReal iterativeResidual(const ADReal & x);
  virtual ADReal iterativeResidualDerivative(const ADReal & x);
  virtual ADReal newtonRoot();
  virtual ADReal brentRoot(const ADReal & x1, const ADReal x2);
  virtual ADReal safeNewtonRoot(const ADReal & x1, const ADReal x2);

  const bool _coupled_temp;
  const ADVariableValue & _temp;

  // Creep parameters
  const bool _has_diffusion_creep;
  const std::vector<Real> _A_diffusion;
  const std::vector<Real> _E_diffusion;
  const std::vector<Real> _V_diffusion;
  const bool _has_dislocation_creep;
  const std::vector<Real> _A_dislocation;
  const std::vector<Real> _n_dislocation;
  const std::vector<Real> _E_dislocation;
  const std::vector<Real> _V_dislocation;
  const Real _gas_constant;
  const std::vector<Real> _eta_min;
  const std::vector<Real> _eta_max;
  const unsigned int _viscous_update;
  const Real _tol;
  const unsigned int _itmax;

  // Creep properties
  ADMaterialProperty(Real) & _eta_eff;
  ADMaterialProperty(RankTwoTensor) & _viscous_strain_incr;

  // Creep effective parameters
  ADReal _A_diff;
  ADReal _E_diff;
  ADReal _V_diff;

  ADReal _A_disl;
  ADReal _n_disl;
  ADReal _E_disl;
  ADReal _V_disl;

  // For iterative update
  ADReal _tau_II_tr;
  ADReal _G;
  

  usingMaterialBaseMembers;
};