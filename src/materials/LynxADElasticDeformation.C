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

#include "LynxADElasticDeformation.h"
#include "LynxADCreepModel.h"
#include "LynxADPlasticModel.h"

registerADMooseObject("LynxApp", LynxADElasticDeformation);

defineADValidParams(
    LynxADElasticDeformation,
    LynxADDeformationBase,
    params.addClassDescription("Class calculating strain and stress for an elastic rheology.");
    // Elastic moduli parameters
    params.addRequiredRangeCheckedParam<std::vector<Real>>("bulk_modulus",
                                                   "bulk_modulus >= 0.0",
                                                   "The drained bulk modulus of the material.");
    params.addRequiredRangeCheckedParam<std::vector<Real>>("shear_modulus",
                                                   "shear_modulus >= 0.0",
                                                   "The shear modulus of the material.");
    // Creep model
    params.addParam<MaterialName>(
        "creep_model",
        "The material object to use for the creep model for a visco-elastic rheology.");
    // Plastic model
    params.addParam<MaterialName>(
        "plastic_model",
        "The material object to use for the plastic model for a elasto-plastic rheology."););

template <ComputeStage compute_stage>
LynxADElasticDeformation<compute_stage>::LynxADElasticDeformation(
    const InputParameters & parameters)
  : LynxADDeformationBase<compute_stage>(parameters),
    // Elastic moduli parameters
    _bulk_modulus(this->getLynxParam("bulk_modulus")),
    _shear_modulus(this->getLynxParam("shear_modulus")),
    // Creep and platic models
    _has_creep(isParamValid("creep_model")),
    _has_plastic(isParamValid("plastic_model")),
    // Elastic properties
    _plith_old(isCoupled("lithostatic_pressure") ? coupledValueOld("lithostatic_pressure") : _zero),
    _elastic_strain_incr(declareADProperty<RankTwoTensor>("elastic_train_increment")),
    _K(declareADProperty<Real>("bulk_modulus")),
    _G(declareADProperty<Real>("shear_modulus")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>("stress")),
    _viscous_strain_incr(_has_creep ? &getADMaterialProperty<RankTwoTensor>("viscous_strain_increment") : nullptr),
    _plastic_strain_incr(_has_plastic ? &getADMaterialProperty<RankTwoTensor>("plastic_strain_increment") : nullptr)
{
}

template <ComputeStage compute_stage>
void
LynxADElasticDeformation<compute_stage>::initialSetup()
{
  LynxADDeformationBase<compute_stage>::initialSetup();
  
  if (_has_creep)
  {
    MaterialName creep_model = getParam<MaterialName>("creep_model");

    LynxADCreepModel<compute_stage> * creep_r = dynamic_cast<LynxADCreepModel<compute_stage> *>(
        &this->template getMaterialByName<compute_stage>(creep_model));

    _creep_model = creep_r;
  }
  else
    _creep_model = nullptr;

  if (_has_plastic)
  {
    MaterialName plastic_model = getParam<MaterialName>("plastic_model");

    LynxADPlasticModel<compute_stage> * plastic_r = dynamic_cast<LynxADPlasticModel<compute_stage> *>(
        &this->template getMaterialByName<compute_stage>(plastic_model));

    _plastic_model = plastic_r;  
  }
  else
    _plastic_model = nullptr;
}

template <ComputeStage compute_stage>
void
LynxADElasticDeformation<compute_stage>::initQpStatefulProperties()
{
  _stress[_qp].zero();
  _stress[_qp].addIa(-_plith[_qp]);
}

template <ComputeStage compute_stage>
void
LynxADElasticDeformation<compute_stage>::initializeQpDeformation()
{
  // Initialize elastic properties
  _elastic_strain_incr[_qp] = _strain_increment[_qp];
  _K[_qp] = this->averageProperty(_bulk_modulus);
  _G[_qp] = this->averageProperty(_shear_modulus);
}

template <ComputeStage compute_stage>
void
LynxADElasticDeformation<compute_stage>::computeQpStress()
{
  // Update the volumetric part of the deformation
  ADReal pressure = volumetricDeformation();

  // Update the deviatoric part of the deformation
  ADRankTwoTensor stress_dev = deviatoricDeformation(pressure);

  // Creep update
  if (_has_creep)
  {
    _creep_model->setQp(_qp);
    _creep_model->creepUpdate(stress_dev, pressure, _G[_qp], _elastic_strain_incr[_qp]);
  }

  // Plastic update
  if (_has_plastic)
  {
    _plastic_model->setQp(_qp);
    _plastic_model->plasticUpdate(stress_dev, pressure, _G[_qp], _K[_qp], _elastic_strain_incr[_qp]);
  }

  // Form the total stress tensor
  LynxADDeformationBase<compute_stage>::reformStressTensor(pressure, stress_dev);
}

template <ComputeStage compute_stage>
ADReal
LynxADElasticDeformation<compute_stage>::volumetricDeformation()
{
  ADReal pressure = -_stress_old[_qp].trace() / 3.0;

  pressure -= _K[_qp] * _elastic_strain_incr[_qp].trace();

  pressure += _plith[_qp] - _plith_old[_qp];

  return pressure;
}

template <ComputeStage compute_stage>
ADRankTwoTensor
LynxADElasticDeformation<compute_stage>::deviatoricDeformation(const ADReal & /*pressure*/)
{
  ADRankTwoTensor stress_dev = this->spinRotation(_stress_old[_qp].deviatoric());

  stress_dev += 2.0 * _G[_qp] * _elastic_strain_incr[_qp].deviatoric();

  return stress_dev; 
}

template <ComputeStage compute_stage>
void
LynxADElasticDeformation<compute_stage>::computeQpThermalSources()
{
  ADRankTwoTensor inelastic_strain_incr = ADRankTwoTensor();
  if (_has_creep)
    inelastic_strain_incr += (*_viscous_strain_incr)[_qp];
  if (_has_plastic)
    inelastic_strain_incr += (*_plastic_strain_incr)[_qp];
  // if (_coupled_temp && !_coupled_temp_aux) 
  //   inelastic_strain_incr.addIa((*_thermal_exp)[_qp] * _temp_dot[_qp] * _dt / 3.0);
  // else if (_coupled_temp_aux && !_coupled_temp)
  //   inelastic_strain_incr.addIa((*_thermal_exp)[_qp] * _temp_dot_aux[_qp] * _dt / 3.0);

  _inelastic_heat[_qp] = _stress[_qp].doubleContraction(inelastic_strain_incr) / _dt;
}