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

#ifndef LYNXMATERIALBASE_H
#define LYNXMATERIALBASE_H

#include "Material.h"
#include "DerivativeMaterialInterface.h"

class LynxMaterialBase;

template <>
InputParameters validParams<LynxMaterialBase>();

class LynxMaterialBase : public DerivativeMaterialInterface<Material>
{
public:
  LynxMaterialBase(const InputParameters & parameters);
  virtual ~LynxMaterialBase() {}
  static MooseEnum averageType();
  template <typename T>
  const std::vector<T> & getLynxParam(const std::string & name) const;

protected:
  virtual Real averageProperty(const std::vector<Real> & properties);
  virtual Real arithmetic_average(const std::vector<Real> & properties);
  virtual Real harmonic_average(const std::vector<Real> & properties);
  virtual Real max_average(const std::vector<Real> & properties);

  bool _has_compositional_phases;
  unsigned int _n_composition;
  const MooseEnum _average_type;
  std::vector<const VariableValue *> _compositional_phases;
};

#endif // LYNXSTRAINBASE_H
