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

#include "LynxADPressureLoad.h"

registerMooseObject("LynxApp", LynxADPressureLoad);

InputParameters
LynxADPressureLoad::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription(
      "Kernel calculating the lithostatic pressure based on density distribution.");
  return params;
}

LynxADPressureLoad::LynxADPressureLoad(const InputParameters & parameters)
  : ADKernel(parameters),
    _bulk_density(getADMaterialProperty<Real>("reference_bulk_density")),
    _gravity(getMaterialProperty<RealVectorValue>("gravity_vector"))
{
}

ADReal
LynxADPressureLoad::computeQpResidual()
{
  return (_grad_u[_qp] - _bulk_density[_qp] * _gravity[_qp]) * _grad_test[_i][_qp];
}