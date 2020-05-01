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

#include "LynxADVolStrainAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LynxApp", LynxADVolStrainAux);

InputParameters
LynxADVolStrainAux::validParams()
{
  InputParameters params = LynxADStrainAuxBase::validParams();
  params.addClassDescription(
      "Access the volumetric part of the strain (total, inelastic, creep or plastic) tensor.");
  return params;
}

LynxADVolStrainAux::LynxADVolStrainAux(const InputParameters & parameters)
  : LynxADStrainAuxBase(parameters)
{
}

Real
LynxADVolStrainAux::computeValue()
{
  return _u_old[_qp] + MetaPhysicL::raw_value((*_strain_incr)[_qp].trace());
}
