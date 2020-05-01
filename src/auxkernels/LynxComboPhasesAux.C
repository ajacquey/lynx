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

#include "LynxComboPhasesAux.h"

registerMooseObject("LynxApp", LynxComboPhasesAux);

InputParameters
LynxComboPhasesAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
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