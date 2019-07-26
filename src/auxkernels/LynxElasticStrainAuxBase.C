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

#include "LynxElasticStrainAuxBase.h"

template <>
InputParameters
validParams<LynxElasticStrainAuxBase>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription(
      "Base class to access the elastic strain tensor.");
  return params;
}

LynxElasticStrainAuxBase::LynxElasticStrainAuxBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<AuxKernel>(parameters),
    _elastic_strain(getDefaultMaterialProperty<RankTwoTensor>("elastic_strain"))
{
}
