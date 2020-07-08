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

#include "LynxADStrainAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LynxApp", LynxADStrainAux);

InputParameters
LynxADStrainAux::validParams()
{
  InputParameters params = LynxADStrainAuxBase::validParams();
  params.addClassDescription(
      "Access a component of the strain (total, inelastic or plastic) tensor.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_i",
      "index_i >= 0 & index_i <= 2",
      "The index i of ij for the stress tensor (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_j",
      "index_j >= 0 & index_j <= 2",
      "The index j of ij for the stress tensor (0, 1, 2)");
  return params;
}

LynxADStrainAux::LynxADStrainAux(const InputParameters & parameters)
  : LynxADStrainAuxBase(parameters),
    _i(getParam<unsigned int>("index_i")),
    _j(getParam<unsigned int>("index_j"))
{
}

Real
LynxADStrainAux::computeValue()
{
  Real strain_incr = MetaPhysicL::raw_value((*_strain_incr)[_qp](_i, _j));
  return (_is_transient) ? _u_old[_qp] + strain_incr : strain_incr;
}
