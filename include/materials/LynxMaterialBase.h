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