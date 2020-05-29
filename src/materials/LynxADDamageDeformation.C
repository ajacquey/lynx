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

#include "LynxADDamageDeformation.h"
#include "ElasticityTensorTools.h"
#include "LynxADCreepModel.h"
#include "LynxADPlasticModel.h"

registerMooseObject("LynxApp", LynxADDamageDeformation);

InputParameters
LynxADDamageDeformation::validParams()
{
  InputParameters params = LynxADElasticDeformation::validParams();
  params.addClassDescription("Class calculating the deformation of a material following a "
                             "damage viscoelasto-(visco)plastic rheology.");
  params.addRequiredParam<MaterialName>(
      "damage_model",
      "The material object to use for the damage rheology.");
  return params;
}

LynxADDamageDeformation::LynxADDamageDeformation(const InputParameters & parameters)
  : LynxADElasticDeformation(parameters),
    // Damage model
    // _damage_model(&getMaterialByName("damage_model")),
    // Strain properties
    _elastic_strain(declareADProperty<RankTwoTensor>("elastic_strain")),
    _elastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("elastic_strain"))
{
  LynxADDamageModelBase * damage_r =
        dynamic_cast<LynxADDamageModelBase *>(&this->getMaterialByName("damage_model"));

    _damage_model = damage_r;
}

void
LynxADDamageDeformation::initializeQpDeformation()
{
  // Initialize elastic properties
  _elastic_strain_incr[_qp] = _strain_increment[_qp];

  _K[_qp] = averageProperty(_bulk_modulus);
  _G[_qp] = averageProperty(_shear_modulus);
  _Cijkl = ElasticityTensorTools::elasticityTensorKandG(_K[_qp], _G[_qp]);
  
  _elastic_strain[_qp] = spinRotation(_elastic_strain_old[_qp]) + _elastic_strain_incr[_qp];
}

void
LynxADDamageDeformation::computeQpStress()
{
  ADRankTwoTensor stress_old = spinRotation(_stress_old[_qp]);
  ADRankTwoTensor elastic_strain_old = spinRotation(_elastic_strain_old[_qp]);
  
  // Damage elastic guess
  _damage_model->setQp(_qp);
  _damage_model->elasticGuess(_stress[_qp], stress_old, _Cijkl, elastic_strain_old, _elastic_strain_incr[_qp]);

  _stress[_qp].addIa(_plith_old[_qp] - _plith[_qp]);

  // Damage - plastic correction

  _damage_model->damageUpdate(_stress[_qp], elastic_strain_old);
}

// void
// LynxADDamageDeformation::computeQpThermalSources()
// {
//   LynxADDeformationBase::computeQpThermalSources();

//   if (_has_plasticity)
//   {
//     Real I2 = Utility::pow<2>(_elastic_strain[_qp].L2norm());
//     Real xi = strainRatio(_elastic_strain[_qp]);
//     Real dPsi_dalpha = _damage_plasticity->_gamma * I2 * (xi - _damage_plasticity->_xi0);
//     _damage_heat[_qp] = dPsi_dalpha * _damage_rate[_qp];
//   }
//   else
//     _damage_heat[_qp] = 0.0;
// }
