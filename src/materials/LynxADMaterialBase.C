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

#include "LynxADMaterialBase.h"

InputParameters
LynxADMaterialBase::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("Base class for compositional phases.");
  params.addCoupledVar("compositional_phases", "The compositional phases.");
  MooseEnum avg_type("arithmetic=0 harmonic=1 max=2", "arithmetic");
  params.addParam<MooseEnum>(
      "average_type",
      avg_type,
      "The type of averaging rule for compositional field dependent properties.");
  return params;
}

LynxADMaterialBase::LynxADMaterialBase(const InputParameters & parameters)
  : ADMaterial(parameters),
    _has_compositional_phases(isCoupled("compositional_phases")),
    _n_composition(_has_compositional_phases ? coupledComponents("compositional_phases") : 1),
    _average_type(getParam<MooseEnum>("average_type")),
    _compositional_phases(_n_composition)
{
  if (_has_compositional_phases)
    for (unsigned i = 0; i < _n_composition; ++i)
      _compositional_phases[i] = &adCoupledValue("compositional_phases", i);
}

template <typename T>
const std::vector<T> &
LynxADMaterialBase::getLynxParam(const std::string & name) const
{
  const std::vector<T> & prop =
      InputParameters::getParamHelper(name, _pars, static_cast<std::vector<T> *>(0));

  if (prop.size() != _n_composition)
    mooseError("Size of vector \"", name, "\" must match the size of compositional phases!");
  else
    return prop;
}

ADReal
LynxADMaterialBase::averageProperty(const std::vector<Real> & properties)
{
  if (_n_composition == 1)
    return properties[0];
  ADReal average_property = 0.0;
  switch (_average_type)
  {
    case 0: // ARITHMETIC
      average_property = arithmetic_average(properties);
      break;
    case 1: // HARMONIC
      average_property = harmonic_average(properties);
      break;
    case 2: // MAX
      average_property = max_average(properties);
      break;
  }

  return average_property;
}

ADReal
LynxADMaterialBase::arithmetic_average(const std::vector<Real> & properties)
{
  ADReal phase = 0;
  ADReal sum = 0;
  ADReal value = 0;
  for (unsigned i = 0; i < _n_composition; ++i)
  {
    phase = std::min(std::max(0.0, (*_compositional_phases[i])[_qp]), 1.0);
    sum += phase;
    value += phase * properties[i];
  }
  value /= sum;

  return value;
}

ADReal
LynxADMaterialBase::harmonic_average(const std::vector<Real> & properties)
{
  ADReal phase = 0;
  ADReal sum = 0;
  ADReal value = 0;
  for (unsigned i = 0; i < _n_composition; ++i)
  {
    phase = std::min(std::max(0.0, (*_compositional_phases[i])[_qp]), 1.0);
    sum += phase;
    value += (properties[i] != 0) ? phase / properties[i] : 0.0;
  }
  value /= sum;

  return (value != 0) ? 1.0 / value : 0;
}

ADReal
LynxADMaterialBase::max_average(const std::vector<Real> & properties)
{
  ADReal max_phase = (*_compositional_phases[0])[_qp];
  ADReal value = properties[0];
  for (unsigned i = 1; i < _n_composition; ++i)
    if ((*_compositional_phases[i])[_qp] > max_phase)
    {
      max_phase = (*_compositional_phases[i])[_qp];
      value = properties[i];
    }

  return value;
}