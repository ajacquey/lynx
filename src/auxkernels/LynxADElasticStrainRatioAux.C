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

#include "LynxADElasticStrainRatioAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LynxApp", LynxADElasticStrainRatioAux);

InputParameters
LynxADElasticStrainRatioAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the elastic strain ratio.");
  return params;
}

LynxADElasticStrainRatioAux::LynxADElasticStrainRatioAux(const InputParameters & parameters)
  : AuxKernel(parameters), _elastic_strain(getADMaterialProperty<RankTwoTensor>("elastic_strain"))
{
}

Real
LynxADElasticStrainRatioAux::computeValue()
{
  Real vol_strain = MetaPhysicL::raw_value(_elastic_strain[_qp].trace());
  Real e_norm = MetaPhysicL::raw_value(_elastic_strain[_qp].L2norm());
  return (e_norm != 0.0) ? vol_strain / e_norm : -std::sqrt(3.0);
}