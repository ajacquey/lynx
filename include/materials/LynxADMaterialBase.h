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

#include "ADMaterial.h"

#define usingMaterialBaseMembers                                                                   \
  usingMaterialMembers;                                                                            \
  using LynxADMaterialBase<compute_stage>::_has_compositional_phases;                              \
  using LynxADMaterialBase<compute_stage>::_n_composition;                                         \
  using LynxADMaterialBase<compute_stage>::_average_type;                                          \
  using LynxADMaterialBase<compute_stage>::_compositional_phases

template <ComputeStage>
class LynxADMaterialBase;

declareADValidParams(LynxADMaterialBase);

template <ComputeStage compute_stage>
class LynxADMaterialBase : public ADMaterial<compute_stage>
{
public:
  LynxADMaterialBase(const InputParameters & parameters);
  const std::vector<Real> & getLynxParam(const std::string & name) const;

protected:
  virtual ADReal averageProperty(const std::vector<Real> & properties);
  virtual ADReal arithmetic_average(const std::vector<Real> & properties);
  virtual ADReal harmonic_average(const std::vector<Real> & properties);
  virtual ADReal max_average(const std::vector<Real> & properties);

  const bool _has_compositional_phases;
  const unsigned int _n_composition;
  const unsigned int  _average_type;
  std::vector<const ADVariableValue *> _compositional_phases;

  usingMaterialMembers;
};