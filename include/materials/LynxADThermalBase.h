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

#define usingThermalBaseMembers                                                                    \
  usingMaterialBaseMembers;                                                                        \
  using LynxADThermalBase<compute_stage>::_c_f;                                                    \
  using LynxADThermalBase<compute_stage>::_c_s;                                                    \
  using LynxADThermalBase<compute_stage>::_lambda_f;                                               \
  using LynxADThermalBase<compute_stage>::_lambda_s;                                               \
  using LynxADThermalBase<compute_stage>::_beta_f;                                                 \
  using LynxADThermalBase<compute_stage>::_beta_s

template <ComputeStage>
class LynxADThermalBase;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;

declareADValidParams(LynxADThermalBase);

template <ComputeStage compute_stage>
class LynxADThermalBase : public LynxADMaterialBase<compute_stage>
{
public:
  LynxADThermalBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void computeQpSpecificHeat();
  virtual void computeQpThermalDiff();
  virtual void computeQpThermalExpansion();
  virtual void computeQpThermalSource();
  virtual void computeQpHeatCap() = 0;
  virtual void computeQpThermalCond() = 0;
  virtual void computeQpThermalExp() = 0;
  virtual ADReal computeMixtureProperty(const ADReal fluid_prop, const ADReal solid_prop);

  // const bool _coupled_porosity;
  const ADVariableValue & _porosity;
  const std::vector<Real> _heat_source;
  const bool _coupled_dens;

  const ADMaterialProperty(Real) * _rho_f;
  const ADMaterialProperty(Real) * _rho_s;
  ADMaterialProperty(Real) & _thermal_diff;
  ADMaterialProperty(Real) & _rhoC_b;
  ADMaterialProperty(Real) & _rhoC_f;
  ADMaterialProperty(Real) & _thermal_exp;
  ADMaterialProperty(Real) & _radiogenic_heat;

  std::vector<ADReal> _c_f;
  std::vector<ADReal> _c_s;
  std::vector<ADReal> _lambda_f;
  std::vector<ADReal> _lambda_s;
  std::vector<ADReal> _beta_f;
  std::vector<ADReal> _beta_s;

  usingMaterialBaseMembers;
};