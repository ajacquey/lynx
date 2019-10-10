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

#include "LynxADDeformationBase.h"

template <ComputeStage>
class LynxADStokeDeformation;

declareADValidParams(LynxADStokeDeformation);

template <ComputeStage compute_stage>
class LynxADStokeDeformation : public LynxADDeformationBase<compute_stage>
{
public:
  LynxADStokeDeformation(const InputParameters & parameters);
protected:
  virtual void initializeQpDeformation() override;
  virtual ADReal volumetricDeformation() override;
  virtual ADRankTwoTensor deviatoricDeformation(const ADReal & pressure) override;
  virtual ADReal computeQpEffectiveViscosity(const ADReal & pressure);
  virtual ADReal computeQpOneOnDiffViscosity(const ADReal A);
  virtual ADReal computeQpOneOnDislViscosity(const ADReal A, const ADReal n, const ADReal eII);
  virtual void computeQpThermalSources() override;

  const bool _coupled_pdyn;
  const ADVariableValue * _pdyn;
  const ADVariableValue * _temp;

  // Stoke parameters
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
  const bool _has_background_strain_rate;
  const bool _has_initial_viscosity;
  const std::vector<Real> _initial_viscosity;
  const Real _background_strain_rate;
  const std::vector<Real> _eta_min;
  const std::vector<Real> _eta_max;

  // Elastic properties
  ADMaterialProperty(Real) & _eta_eff;

  ADReal _A_diff;
  ADReal _E_diff;
  ADReal _V_diff;

  ADReal _A_disl;
  ADReal _n_disl;
  ADReal _E_disl;
  ADReal _V_disl;

  usingDeformationBaseMembers;
};