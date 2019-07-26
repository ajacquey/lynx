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

#ifndef LYNXELASTICSTRAINAUXBASE_H
#define LYNXELASTICSTRAINAUXBASE_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "DerivativeMaterialInterface.h"

class LynxElasticStrainAuxBase;

template <>
InputParameters validParams<LynxElasticStrainAuxBase>();

class LynxElasticStrainAuxBase : public DerivativeMaterialInterface<AuxKernel>
{
public:
  LynxElasticStrainAuxBase(const InputParameters & parameters);

  const MaterialProperty<RankTwoTensor> & _elastic_strain;
};

#endif // LYNXELASTICSTRAINAUXBASE_H
