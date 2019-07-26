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

#include "LynxElasticVolStrainAux.h"

registerMooseObject("LynxApp", LynxElasticVolStrainAux);

template <>
InputParameters
validParams<LynxElasticVolStrainAux>()
{
  InputParameters params = validParams<LynxElasticStrainAuxBase>();
  params.addClassDescription(
      "Access the volumetric elastic strain.");
  return params;
}

LynxElasticVolStrainAux::LynxElasticVolStrainAux(const InputParameters & parameters)
  : LynxElasticStrainAuxBase(parameters)
{
}

Real
LynxElasticVolStrainAux::computeValue()
{
  return _elastic_strain[_qp].trace();
}
