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

#include "LynxStrainAuxBase.h"

template <>
InputParameters
validParams<LynxStrainAuxBase>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Base class for outputting strain values.");
  params.addParam<MooseEnum>("strain_type",
                             LynxStrainAuxBase::strainType() = "total",
                             "The type of the strain tensor to output.");
  return params;
}

LynxStrainAuxBase::LynxStrainAuxBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<AuxKernel>(parameters),
    _strain_type(getParam<MooseEnum>("strain_type"))
{
  switch (_strain_type)
  {
    case 1:
      _strain_name = "strain_increment";
      break;
    case 2:
      _strain_name = "inelastic_strain_increment";
      break;
    case 3:
      _strain_name = "elastic_strain_increment";
      break;
    case 4:
      _strain_name = "viscous_strain_increment";
      break;
    case 5:
      _strain_name = "damage_strain_increment";
      break;
    case 6:
      _strain_name = "plastic_strain_increment";
      break;
    case 7:
      _strain_name = "thermal_strain_increment";
      break;
    default:
      mooseError("LynxStrainAuxBase: unknown strain type!");
  }
  _strain_incr = &getDefaultMaterialProperty<RankTwoTensor>(_strain_name);
}

MooseEnum
LynxStrainAuxBase::strainType()
{
  return MooseEnum("total=1 inelastic=2 elastic=3 viscous=4 damage=5 plastic=6 thermal=7");
}
