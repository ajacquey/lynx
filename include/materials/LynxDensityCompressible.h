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
#ifndef LYNXDENSITYCOMPRESSIBLE_H
#define LYNXDENSITYCOMPRESSIBLE_H

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

#endif // LYNXDENSITYCOMPRESSIBLE_H
