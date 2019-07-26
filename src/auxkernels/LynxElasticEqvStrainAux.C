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

#include "LynxElasticEqvStrainAux.h"

registerMooseObject("LynxApp", LynxElasticEqvStrainAux);

template <>
InputParameters
validParams<LynxElasticEqvStrainAux>()
{
  InputParameters params = validParams<LynxElasticStrainAuxBase>();
  params.addClassDescription(
      "Access the volumetric elastic strain.");
  return params;
}

LynxElasticEqvStrainAux::LynxElasticEqvStrainAux(const InputParameters & parameters)
  : LynxElasticStrainAuxBase(parameters)
{
}

Real
LynxElasticEqvStrainAux::computeValue()
{
  return std::sqrt(2.0 / 3.0) * _elastic_strain[_qp].deviatoric().L2norm();
}
