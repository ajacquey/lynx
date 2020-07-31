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

#include "LynxADElasticDeformation.h"
#include "LynxADDamageModelBase.h"

class LynxADDamageDeformation : public LynxADElasticDeformation
{
public:
  static InputParameters validParams();
  LynxADDamageDeformation(const InputParameters & parameters);
  void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void initializeQpDeformation() override;
  virtual void computeQpStress() override;

  // Damage rheology
  LynxADDamageModelBase * _damage_model;

  // Strain properties
  ADMaterialProperty<RankTwoTensor> & _elastic_strain;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;

  // Elasticity tensor
  ADRankFourTensor _Cijkl;
};