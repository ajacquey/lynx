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

#define usingDeformationBaseMembers                                                                \
  usingMaterialBaseMembers;                                                                        \
  using LynxADDeformationBase<compute_stage>::_plith;                                              \
  using LynxADDeformationBase<compute_stage>::_coupled_temp;                                       \
  using LynxADDeformationBase<compute_stage>::_temp_dot;                                           \
  using LynxADDeformationBase<compute_stage>::_coupled_temp_aux;                                   \
  using LynxADDeformationBase<compute_stage>::_temp_dot_aux;                                       \
  using LynxADDeformationBase<compute_stage>::_strain_increment;                                   \
  using LynxADDeformationBase<compute_stage>::_spin_increment;                                     \
  using LynxADDeformationBase<compute_stage>::_thermal_exp;                                        \
  using LynxADDeformationBase<compute_stage>::_stress;                                             \
  using LynxADDeformationBase<compute_stage>::_inelastic_heat;                                     \
  using LynxADDeformationBase<compute_stage>::spinRotation

template <ComputeStage>
class LynxADDeformationBase;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;

declareADValidParams(LynxADDeformationBase);

template <ComputeStage compute_stage>
class LynxADDeformationBase : public LynxADMaterialBase<compute_stage>
{
public:
  LynxADDeformationBase(const InputParameters & parameters);
  void initialSetup() override;
  void displacementIntegrityCheck();

protected:
  virtual void computeProperties() override;
  virtual void computeQpProperties() override;
  virtual void computeStrainIncrement();
  virtual void calculateSmallStrain(const ADRankTwoTensor & grad_tensor,
                                    const RankTwoTensor & grad_tensor_old);
  virtual void calculateFiniteStrain(const ADRankTwoTensor & grad_tensor,
                                     const RankTwoTensor & grad_tensor_old);
  virtual void computeQpDeformation();
  virtual void initializeQpDeformation();
  virtual void computeQpStress();
  virtual ADReal volumetricDeformation() = 0;
  virtual ADRankTwoTensor deviatoricDeformation(const ADReal & pressure) = 0;
  virtual void reformStressTensor(const ADReal & pressure, const ADRankTwoTensor & stress_dev);
  virtual ADRankTwoTensor spinRotation(const ADRankTwoTensor & tensor);
  virtual void computeQpThermalSources() = 0;

  // Coupled variables
  const unsigned int _ndisp;
  std::vector<const ADVariableGradient *> _grad_disp;
  std::vector<const VariableGradient *> _grad_disp_old;
  const ADVariableValue & _plith;
  const bool _coupled_temp;
  const ADVariableValue & _temp_dot;
  const bool _coupled_temp_aux;
  const ADVariableValue & _temp_dot_aux;

  // Strain parameters
  const unsigned int _strain_model;
  const bool _vol_locking_correction;
  const Real & _current_elem_volume;

  // Strain properties
  ADMaterialProperty(RankTwoTensor) & _strain_increment;
  ADMaterialProperty(RankTwoTensor) & _spin_increment;
  const ADMaterialProperty(Real) * _thermal_exp;

  // Stress properties
  ADMaterialProperty(RankTwoTensor) & _stress;
  ADMaterialProperty(Real) & _inelastic_heat;

  usingMaterialBaseMembers;
};