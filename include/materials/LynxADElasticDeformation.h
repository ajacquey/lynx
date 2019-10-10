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
class LynxADElasticDeformation;
template <ComputeStage>
class LynxADCreepModel;
template <ComputeStage>
class LynxADPlasticModel;

declareADValidParams(LynxADElasticDeformation);

template <ComputeStage compute_stage>
class LynxADElasticDeformation : public LynxADDeformationBase<compute_stage>
{
public:
  LynxADElasticDeformation(const InputParameters & parameters);
  void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void initializeQpDeformation() override;
  virtual void computeQpStress() override;
  virtual ADReal volumetricDeformation() override;
  virtual ADRankTwoTensor deviatoricDeformation(const ADReal & pressure) override;
  virtual void computeQpThermalSources() override;

  // Elastic parameters
  const std::vector<Real> _bulk_modulus;
  const std::vector<Real> _shear_modulus;

  // Creep and plastic models
  const bool _has_creep;
  const bool _has_plastic;

  // Elastic properties
  const VariableValue & _plith_old;
  ADMaterialProperty(RankTwoTensor) & _elastic_strain_incr;
  ADMaterialProperty(Real) & _K;
  ADMaterialProperty(Real) & _G;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  // Creep Model
  LynxADCreepModel<compute_stage> * _creep_model;
  // Plastic Model
  LynxADPlasticModel<compute_stage> * _plastic_model;

  usingDeformationBaseMembers;
};