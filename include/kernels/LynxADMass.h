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

#include "ADKernel.h"

class LynxADMass : public ADKernel
{
public:
  static InputParameters validParams();
  LynxADMass(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _penalty;
  const unsigned int _penalty_type;
  const ADMaterialProperty<RankTwoTensor> & _strain_increment;
};