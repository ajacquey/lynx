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

#include "LynxDensityBase.h"

class LynxDensityCompressible;

template <>
InputParameters validParams<LynxDensityCompressible>();

class LynxDensityCompressible : public LynxDensityBase
{
public:
  LynxDensityCompressible(const InputParameters & parameters);

protected:
  void computeQpProperties();

  const bool _coupled_temp;
  const VariableValue & _temp;
  const VariableValue & _reference_temperature;
  const std::vector<Real> _beta_solid;
  // const bool _temperature_from_multiapp;
  MaterialProperty<Real> & _drho_dtemp;
  MaterialProperty<Real> & _drho_dev;
  MaterialProperty<Real> & _dinvrho_dtemp;
  MaterialProperty<Real> & _dinvrho_dev;

  const bool _coupled_plith;
  const VariableValue & _plith;
  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<Real> & _K;
};