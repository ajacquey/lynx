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

#include "LynxADMaterialBase.h"

#define usingDensityBaseMembers                                                                    \
  usingMaterialBaseMembers;                                                                        \
  using LynxADDensityBase<compute_stage>::_porosity;                                               \
  using LynxADDensityBase<compute_stage>::_fluid_density;                                          \
  using LynxADDensityBase<compute_stage>::_solid_density;                                          \
  using LynxADDensityBase<compute_stage>::_rho_f;                                                  \
  using LynxADDensityBase<compute_stage>::_rho_s;                                                  \
  using LynxADDensityBase<compute_stage>::_rho_b;                                                  \
  using LynxADDensityBase<compute_stage>::_reference_rho_b;                                        \
  using LynxADDensityBase<compute_stage>::computeQpGravity

template <ComputeStage>
class LynxADDensityBase;

declareADValidParams(LynxADDensityBase);

template <ComputeStage compute_stage>
class LynxADDensityBase : public LynxADMaterialBase<compute_stage>
{
public:
  LynxADDensityBase(const InputParameters & parameters);

protected:
  virtual void computeQpGravity();

  const ADVariableValue & _porosity;

  bool _has_gravity;
  Real _g;
  const std::vector<Real> _fluid_density;
  const std::vector<Real> _solid_density;

  ADMaterialProperty(RealVectorValue) & _gravity;
  ADMaterialProperty(Real) & _rho_f;
  ADMaterialProperty(Real) & _rho_s;
  ADMaterialProperty(Real) & _rho_b;
  ADMaterialProperty(Real) & _reference_rho_b;

  usingMaterialBaseMembers;
};