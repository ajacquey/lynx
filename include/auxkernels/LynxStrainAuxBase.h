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

#ifndef LYNXSTRAINAUXBASE_H
#define LYNXSTRAINAUXBASE_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "DerivativeMaterialInterface.h"

class LynxStrainAuxBase;

template <>
InputParameters validParams<LynxStrainAuxBase>();

class LynxStrainAuxBase : public DerivativeMaterialInterface<AuxKernel>
{
public:
  LynxStrainAuxBase(const InputParameters & parameters);
  static MooseEnum strainType();

protected:
  MooseEnum _strain_type;
  std::string _strain_name;
  const MaterialProperty<RankTwoTensor> * _strain_incr;
};

#endif // LYNXSTRAINAUXBASE_H