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

#ifndef LYNXDEFORMATIONBASE_H
#define LYNXDEFORMATIONBASE_H

#include "LynxMaterialBase.h"
#include "LynxRheologyStructures.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class LynxDeformationBase;

template <>
InputParameters validParams<LynxDeformationBase>();

class LynxDeformationBase : public LynxMaterialBase
{
public:
  LynxDeformationBase(const InputParameters & parameters);
  virtual ~LynxDeformationBase() {}
  static MooseEnum strainModel();
  static MooseEnum viscousUpdate();

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeProperties() override;
  virtual void computeQpProperties() override;
  virtual void computeStrainIncrement();
  virtual void calculateSmallStrain(const RankTwoTensor & grad_tensor,
                                    const RankTwoTensor & grad_tensor_old);
  virtual void calculateFiniteStrain(const RankTwoTensor & grad_tensor,
                                     const RankTwoTensor & grad_tensor_old);
  virtual void computeQpDeformation();
  virtual void initializeQpDeformation();
  virtual void elasticModuli();
  virtual void volumetricDeformation(Real & pressure);
  virtual void deviatoricDeformation(const Real & pressure, RankTwoTensor & stress_dev);
  virtual void reformStressTensor(const Real & pressure, const RankTwoTensor & stress_dev);
  virtual void plasticCorrection(Real & pressure, RankTwoTensor & stress_dev) = 0;
  virtual void damageCorrection();
  virtual void tangentOperator();
  virtual void plasticTangentOperator(const RankTwoTensor & flow_direction,
                                      const RankFourTensor & flow_direction_dyad);
  virtual void finiteTangentOperator();
  virtual Real computeStokeEffectiveViscosity(const Real & pressure);
  virtual RankTwoTensor computeStokeEffectiveViscosityDerivative();
  virtual Real computeEffectiveViscosity(const Real & pressure, const Real & eqv_stress);
  virtual Real computeCreepRate(const Real A, const Real n, const Real eqv_stress);
  virtual Real viscousIncrement(const Real & pressure, const Real & eqv_stress);
  virtual RankTwoTensor spinRotation(const RankTwoTensor & tensor);
  virtual void updateSpinTangentModulus();
  virtual void computeCreepProperties(const Real & pressure);
  virtual void updateCreepParameters();
  virtual Real rootBrent(iterative_viscous & viscous_model, const Real x1, const Real x2);
  virtual Real rootNewtonSafe(iterative_viscous & viscous_model, const Real x1, const Real x2);
  void computeQpThermalSources();

  // Coupled variables
  unsigned int _ndisp;
  std::vector<const VariableGradient *> _grad_disp;
  std::vector<const VariableGradient *> _grad_disp_old;
  bool _coupled_temp;
  const VariableValue & _temp;
  bool _coupled_plith;
  const VariableValue & _plith;
  const VariableValue & _plith_old;
  bool _coupled_pdyn;
  const VariableValue & _pdyn;

  // Strain parameters
  const MooseEnum _strain_model;
  bool _vol_locking_correction;
  std::vector<RankTwoTensor> _deformation_gradient;
  std::vector<RankTwoTensor> _deformation_gradient_old;
  const Real & _current_elem_volume;

  // Elastic moduli parameters
  const bool _has_bulk_modulus;
  const bool _has_shear_modulus;
  const std::vector<Real> _bulk_modulus;
  const std::vector<Real> _shear_modulus;

  // Creep parameters
  const bool _has_diffusion_creep;
  const std::vector<Real> _A_diffusion;
  const std::vector<Real> _E_diffusion;
  const std::vector<Real> _V_diffusion;
  const bool _has_dislocation_creep;
  const std::vector<Real> _A_dislocation;
  const std::vector<Real> _E_dislocation;
  const std::vector<Real> _V_dislocation;
  const std::vector<Real> _n_dislocation;
  const Real _gas_constant;
  const MooseEnum _viscous_update;
  const Real _tol;
  const unsigned int _itmax;
  const bool _has_background_strain_rate;
  bool _has_initial_viscosity;
  const std::vector<Real> _initial_viscosity;
  const Real _background_strain_rate;
  const std::vector<Real> _eta_min;
  const std::vector<Real> _eta_max;

  // Rheology boolean
  const bool _has_elastic;
  const bool _has_viscous;
  bool _has_plasticity;

  // Creep structures
  diffusion_creep * _diffusion_creep;
  dislocation_creep * _dislocation_creep;
  iterative_viscous * _iterative_viscous;

  // RankFourTensor utilities
  RankTwoTensor _identity_two;
  RankFourTensor _volumetric_four;
  RankFourTensor _identity_four;
  RankFourTensor _deviatoric_four;

  // Viscous scalars for tangent modulus
  Real _stress_corr_v;
  Real _dq_dq_tr_v;

  // Plastic scalars for tangent modulus
  Real _stress_corr_p;
  Real _dp_dp_tr_p;
  Real _dp_dq_tr_p;
  Real _dq_dp_tr_p;
  Real _dq_dq_tr_p;

  // Strain properties
  MaterialProperty<RankTwoTensor> & _strain_increment;
  MaterialProperty<RankTwoTensor> & _spin_tensor;
  const MaterialProperty<RankTwoTensor> & _thermal_strain_incr;

  // Viscous properties
  MaterialProperty<RankTwoTensor> & _viscous_strain_incr;
  MaterialProperty<Real> & _eta_eff;

  // Plastic properties
  MaterialProperty<RankTwoTensor> & _plastic_strain_incr;
  MaterialProperty<Real> & _plastic_yield_function;

  // Stress properties
  MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  MaterialProperty<Real> & _K;
  MaterialProperty<Real> & _G;
  MaterialProperty<RankFourTensor> & _tangent_modulus;
  MaterialProperty<Real> & _inelastic_heat;
  MaterialProperty<RankTwoTensor> & _dinelastic_heat_dstrain;
};

#endif // LYNXDEFORMATIONBASE_H
