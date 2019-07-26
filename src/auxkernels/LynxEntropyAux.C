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
