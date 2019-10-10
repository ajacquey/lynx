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
class LynxADHydroDarcy;

declareADValidParams(LynxADHydroDarcy);

template <ComputeStage compute_stage>
class LynxADHydroDarcy : public ADKernel<compute_stage>
{
public:
  LynxADHydroDarcy(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty(Real) & _fluid_mobility;
  const bool _coupled_grav;
  const ADMaterialProperty(RealVectorValue) * _gravity;
  const ADMaterialProperty(Real) * _rho_f;

  usingKernelMembers;
};