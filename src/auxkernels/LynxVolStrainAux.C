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

#include "LynxVolStrainAux.h"

registerMooseObject("LynxApp", LynxVolStrainAux);

template <>
InputParameters
validParams<LynxVolStrainAux>()
{
  InputParameters params = validParams<LynxStrainAuxBase>();
  params.addClassDescription(
      "Access the volumetric part of the strain (total, inelastic, creep or plastic) tensor.");
  return params;
}

LynxVolStrainAux::LynxVolStrainAux(const InputParameters & parameters)
  : LynxStrainAuxBase(parameters)
{
}

Real
LynxVolStrainAux::computeValue()
{
  return _u_old[_qp] + (*_strain_incr)[_qp].trace();
}