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

#include "LynxEqvStrainRateAux.h"

registerMooseObject("LynxApp", LynxEqvStrainRateAux);

template <>
InputParameters
validParams<LynxEqvStrainRateAux>()
{
  InputParameters params = validParams<LynxStrainAuxBase>();
  params.addClassDescription("Calculates the equivalent strain rate of the given tensor.");
  return params;
}

LynxEqvStrainRateAux::LynxEqvStrainRateAux(const InputParameters & parameters)
  : LynxStrainAuxBase(parameters)
{
}

Real
LynxEqvStrainRateAux::computeValue()
{
  return std::sqrt(2.0 / 3.0) * (*_strain_incr)[_qp].deviatoric().L2norm() / _dt;
}
