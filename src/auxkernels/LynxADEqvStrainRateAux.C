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

#include "LynxADEqvStrainRateAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LynxApp", LynxADEqvStrainRateAux);

InputParameters
LynxADEqvStrainRateAux::validParams()
{
  InputParameters params = LynxADStrainAuxBase::validParams();
  params.addClassDescription("Calculates the equivalent strain rate of the given tensor.");
  return params;
}

LynxADEqvStrainRateAux::LynxADEqvStrainRateAux(const InputParameters & parameters)
  : LynxADStrainAuxBase(parameters)
{
}

Real
LynxADEqvStrainRateAux::computeValue()
{
  return std::sqrt(2.0 / 3.0) * MetaPhysicL::raw_value((*_strain_incr)[_qp].deviatoric().L2norm()) /
         _dt;
}
