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
