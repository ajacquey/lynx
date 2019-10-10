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

template <ComputeStage>
class LynxADHydroPoroMech;

declareADValidParams(LynxADHydroPoroMech);

template <ComputeStage compute_stage>
class LynxADHydroPoroMech : public ADKernel<compute_stage>
{
public:
  LynxADHydroPoroMech(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty(Real) & _poro_mech;

  usingKernelMembers;
};