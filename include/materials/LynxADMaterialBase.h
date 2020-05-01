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

class LynxADMaterialBase : public ADMaterial
{
public:
  static InputParameters validParams();
  LynxADMaterialBase(const InputParameters & parameters);
  template <typename T>
  const std::vector<T> & getLynxParam(const std::string & name) const;

protected:
  virtual ADReal averageProperty(const std::vector<Real> & properties);
  virtual ADReal arithmetic_average(const std::vector<Real> & properties);
  virtual ADReal harmonic_average(const std::vector<Real> & properties);
  virtual ADReal max_average(const std::vector<Real> & properties);

  const bool _has_compositional_phases;
  const unsigned int _n_composition;
  const unsigned int  _average_type;
  std::vector<const ADVariableValue *> _compositional_phases;
};