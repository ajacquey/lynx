/******************************************************************************/
/*                       LYNX, a MOOSE-based application                      */
/*                                                                            */
/*          Copyright (C) 2017 by Antoine B. Jacquey and Mauro Cacace         */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*    This program is free software: you can redistribute it and/or modify    */
/*    it under the terms of the GNU General Public License as published by    */
/*      the Free Software Foundation, either version 3 of the License, or     */
/*                     (at your option) any later version.                    */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program. If not, see <http://www.gnu.org/licenses/>     */
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