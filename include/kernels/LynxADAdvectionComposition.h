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

#include "LynxADAdvectionBase.h"

class LynxADAdvectionComposition : public LynxADAdvectionBase
{
public:
  static InputParameters validParams();
  LynxADAdvectionComposition(const InputParameters & parameters);

protected:
  virtual ADReal computeArtificialViscosity() override;
  virtual void computeEntropyResidual();
};