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

#include "LynxADHeatConduction.h"

registerMooseObject("LynxApp", LynxADHeatConduction);

InputParameters
LynxADHeatConduction::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Heat conduction kernel.");
  params.addParam<bool>(
      "use_displaced_mesh", true, "Set the displaced mesh flag to true by default.");
  return params;
}

LynxADHeatConduction::LynxADHeatConduction(const InputParameters & parameters)
  : ADKernel(parameters), _thermal_diff(getADMaterialProperty<Real>("thermal_diffusivity"))
{
}

ADReal
LynxADHeatConduction::computeQpResidual()
{
  return _thermal_diff[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
}