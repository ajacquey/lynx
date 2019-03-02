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

#include "LynxEntropyAux.h"

registerMooseObject("LynxApp", LynxEntropyAux);

template <>
InputParameters
validParams<LynxEntropyAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Computes the entropy of the system needed in the definition of the "
                             "artificial viscosity to stabilize the advection equation - either "
                             "for temperature or any compositional field.");
  params.addRequiredCoupledVar("entropy_variable", "The variable that is advected.");
  params.addRequiredParam<PostprocessorName>(
      "pp_max_var",
      "The postprocessor to retrieve the maximum of the advected variable on the whole domain.");
  params.addRequiredParam<PostprocessorName>(
      "pp_min_var",
      "The postprocessor to retrieve the minimum of the advected variable on the whole domain.");
  return params;
}

LynxEntropyAux::LynxEntropyAux(const InputParameters & parameters)
  : DerivativeMaterialInterface<AuxKernel>(parameters),
    _var_old(_is_transient ? coupledValueOld("entropy_variable") : _zero),
    _var_older(_is_transient ? coupledValueOlder("entropy_variable") : _zero),
    _pp_max_var(getPostprocessorValue("pp_max_var")),
    _pp_min_var(getPostprocessorValue("pp_min_var"))
{
}

Real
LynxEntropyAux::computeValue()
{
  Real var_avg = 0.5 * (_pp_max_var + _pp_min_var);
  Real extra_var = 0.5 * (_var_old[_qp] + _var_older[_qp]);
  return (extra_var - var_avg) * (extra_var - var_avg);
}
