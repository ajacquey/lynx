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

#include "LynxStrainRatioAux.h"

registerMooseObject("LynxApp", LynxStrainRatioAux);

template <>
InputParameters
validParams<LynxStrainRatioAux>()
{
  InputParameters params = validParams<LynxElasticStrainAuxBase>();
  params.addClassDescription(
      "Access the strain ratio (volumetric strain on norm of the elastic strain).");
  return params;
}

LynxStrainRatioAux::LynxStrainRatioAux(const InputParameters & parameters)
  : LynxElasticStrainAuxBase(parameters)
{
}

Real
LynxStrainRatioAux::computeValue()
{
  Real strain_norm = _elastic_strain[_qp].L2norm();

  if (strain_norm != 0.0)
    return _elastic_strain[_qp].trace() / strain_norm;
  else
    return -std::sqrt(3.0);
}
