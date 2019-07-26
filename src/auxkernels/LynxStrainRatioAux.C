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

#include "LynxStrainRatioAux.h"

registerMooseObject("LynxApp", LynxStrainRatioAux);

template <>
InputParameters
validParams<LynxStrainRatioAux>()
{
  InputParameters params = validParams<LynxElasticStrainAuxBase>();
  params.addClassDescription(
      "Access the strain ratio (volumetric strain on norm of the elastic strain).");
  return params;
}

LynxStrainRatioAux::LynxStrainRatioAux(const InputParameters & parameters)
  : LynxElasticStrainAuxBase(parameters)
{
}

Real
LynxStrainRatioAux::computeValue()
{
  Real strain_norm = _elastic_strain[_qp].L2norm();

  if (strain_norm != 0.0)
    return _elastic_strain[_qp].trace() / strain_norm;
  else
    return -std::sqrt(3.0);
}
