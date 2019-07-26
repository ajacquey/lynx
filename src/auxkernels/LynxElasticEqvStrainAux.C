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

#include "LynxElasticEqvStrainAux.h"

registerMooseObject("LynxApp", LynxElasticEqvStrainAux);

template <>
InputParameters
validParams<LynxElasticEqvStrainAux>()
{
  InputParameters params = validParams<LynxElasticStrainAuxBase>();
  params.addClassDescription(
      "Access the volumetric elastic strain.");
  return params;
}

LynxElasticEqvStrainAux::LynxElasticEqvStrainAux(const InputParameters & parameters)
  : LynxElasticStrainAuxBase(parameters)
{
}

Real
LynxElasticEqvStrainAux::computeValue()
{
  return std::sqrt(2.0 / 3.0) * _elastic_strain[_qp].deviatoric().L2norm();
}
