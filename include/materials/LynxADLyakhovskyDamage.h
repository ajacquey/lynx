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

#include "LynxADDamageModelBase.h"

class LynxADLyakhovskyDamage : public LynxADDamageModelBase
{
public:
  static InputParameters validParams();
  LynxADLyakhovskyDamage(const InputParameters & parameters);
  virtual void elasticGuess(ADRankTwoTensor & stress,
                            const ADRankTwoTensor & stress_old,
                            const RankFourTensor & Cijkl,
                            const ADRankTwoTensor & elastic_strain_old,
                            const ADRankTwoTensor & elastic_strain_incr) override;
  virtual void damageUpdate(ADRankTwoTensor & stress,
                            ADRankTwoTensor & elastic_strain_incr) override;

protected:
  virtual ADRankFourTensor damageStiffness(const RankFourTensor & Cijkl);

  virtual void returnMap(ADReal & gamma_s, ADReal & gamma_d);
  virtual void
  residual(const ADReal & gamma_vs, const ADReal & gamma_d, ADReal & ress, ADReal & resd);
  virtual void jacobian(const ADReal & gamma_s,
                        const ADReal & gamma_d,
                        ADReal & jacss,
                        ADReal & jacdd,
                        ADReal & jacsd,
                        ADReal & jacds);
  virtual void
  overStress(const ADReal & gamma_s, const ADReal & gamma_d, ADReal & over_s, ADReal & over_d);
  virtual void overStressDerivS(const ADReal & gamma_s,
                                const ADReal & gamma_d,
                                ADReal & over_s_s,
                                ADReal & over_d_s);
  virtual void overStressDerivD(const ADReal & gamma_s,
                                const ADReal & gamma_d,
                                ADReal & over_s_d,
                                ADReal & over_d_d);
  virtual void updateYieldParameters(const ADReal & gamma_s, const ADReal & gamma_d);
  virtual void updateYieldParametersDerivS(const ADReal & gamma_s, const ADReal & gamma_d);
  virtual void updateYieldParametersDerivD(const ADReal & gamma_s, const ADReal & gamma_d);
  virtual ADReal flowRate();
  virtual ADReal flowDirectionS();
  virtual ADReal flowDirectionD();
  virtual ADReal dFlowRatedS();
  virtual ADReal dFlowRatedD();
  virtual ADReal dFlowDirectionSdS();
  virtual ADReal dFlowDirectionSdD();
  virtual ADReal dFlowDirectionDdS();
  virtual ADReal dFlowDirectionDdD();
  virtual ADReal yieldFunction();
  virtual ADReal dYieldFunctiondDev();
  virtual ADReal dYieldFunctiondPres();
  virtual ADReal dYieldFunctiondDam();
  virtual ADReal dYieldFunctiondMua();
  virtual ADReal d2YieldFunctiondDev2();
  virtual ADReal d2YieldFunctiondDevdPres();
  virtual ADReal d2YieldFunctiondDevdMua();
  virtual ADReal d2YieldFunctiondDamdDev();
  virtual ADReal d2YieldFunctiondDamdPres();
  virtual ADReal d2YieldFunctiondDamdMua();
  virtual ADReal d2YieldFunctiondDam2();
  virtual ADReal strainRatio(const ADRankTwoTensor & elastic_strain);
  virtual ADRankTwoTensor reformPlasticStrainTensor(const ADReal & gamma_s);

  // Damage parameters
  const std::vector<Real> _damage_modulus;
  const std::vector<Real> _friction_angle;
  const std::vector<Real> _cohesion;
  const std::vector<Real> _plastic_viscosity;

  // Damage utilities
  Real _Gam0;
  Real _eta_p;
  Real _K;
  Real _G;
  Real _xi0;
  Real _k;
  ADReal _e_norm;
  ADReal _xi;
  ADRankTwoTensor _e;
  Real _mu0;
  ADRankFourTensor _Dijkl;

  // Return map utilities
  // Trial states
  ADRankTwoTensor _stress_tr;
  ADReal _pressure_tr;
  ADReal _dev_stress_tr;
  ADReal _damage_force_tr;
  // Forces
  ADReal _pressure;
  ADReal _dev_stress;
  ADReal _damage_force;
  ADReal _mua;
  // Forces derivatives
  ADReal _dpressure_dS;
  ADReal _dpressure_dD;
  ADReal _ddev_stress_dS;
  ADReal _ddev_stress_dD;
  ADReal _ddamage_force_dS;
  ADReal _ddamage_force_dD;
  ADReal _dmua_dS;
  ADReal _dmua_dD;
};