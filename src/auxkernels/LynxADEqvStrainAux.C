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

#include "LynxADEqvStrainAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LynxApp", LynxADEqvStrainAux);

InputParameters
LynxADEqvStrainAux::validParams()
{
  InputParameters params = LynxADStrainAuxBase::validParams();
  params.addClassDescription("Calculates the equivalent strain of the given tensor.");
  return params;
}

LynxADEqvStrainAux::LynxADEqvStrainAux(const InputParameters & parameters)
  : LynxADStrainAuxBase(parameters)
{
}

Real
LynxADEqvStrainAux::computeValue()
{
  Real eqv_strain_incr = std::sqrt(2.0 / 3.0) * MetaPhysicL::raw_value((*_strain_incr)[_qp].deviatoric().L2norm());
  return (_is_transient) ? _u_old[_qp] + eqv_strain_incr : eqv_strain_incr;
}
