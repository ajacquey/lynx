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