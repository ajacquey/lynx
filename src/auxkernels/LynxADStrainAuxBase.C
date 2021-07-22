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

#include "LynxADStrainAuxBase.h"

InputParameters
LynxADStrainAuxBase::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Base class for outputting strain values.");
  params.addParam<MooseEnum>("strain_type",
                             LynxADStrainAuxBase::strainType() = "total",
                             "The type of the strain tensor to output.");
  return params;
}

LynxADStrainAuxBase::LynxADStrainAuxBase(const InputParameters & parameters)
  : AuxKernel(parameters),
    _u_old(uOld()),
    _strain_type(getParam<MooseEnum>("strain_type"))
{
  switch (_strain_type)
  {
    case 1:
      _strain_name = "strain_increment";
      break;
    case 2:
      _strain_name = "elastic_strain_increment";
      break;
    case 3:
      _strain_name = "viscous_strain_increment";
      break;
    case 4:
      _strain_name = "plastic_strain_increment";
      break;
    default:
      mooseError("LynxADStrainAuxBase: unknown strain type!");
  }
  _strain_incr = &getADMaterialProperty<RankTwoTensor>(_strain_name);
}

MooseEnum
LynxADStrainAuxBase::strainType()
{
  return MooseEnum("total=1 elastic=2 viscous=3 plastic=4");
}
