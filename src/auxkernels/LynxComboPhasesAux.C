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

#include "LynxComboPhasesAux.h"

registerMooseObject("LynxApp", LynxComboPhasesAux);

template <>
InputParameters
validParams<LynxComboPhasesAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Sample all compositional phases in one aux variable. It gets the "
                             "index of the maximum phase at one quadrature point.");
  params.addCoupledVar("compositional_phases", "The active compositional phases.");
  return params;
}

LynxComboPhasesAux::LynxComboPhasesAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _n_composition(coupledComponents("compositional_phases")),
    _compositional_phases(_n_composition)
{
  for (unsigned i = 0; i < _n_composition; ++i)
    _compositional_phases[i] = &coupledValue("compositional_phases", i);
}

Real
LynxComboPhasesAux::computeValue()
{
  Real phase = (*_compositional_phases[0])[_qp];
  Real value = 0;
  for (unsigned i = 1; i < _n_composition; ++i)
    if ((*_compositional_phases[i])[_qp] > phase)
    {
      phase = (*_compositional_phases[i])[_qp];
      value = i;
    }
  return value;
}
