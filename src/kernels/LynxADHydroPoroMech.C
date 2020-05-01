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

#include "LynxADHydroPoroMech.h"

registerMooseObject("LynxApp", LynxADHydroPoroMech);

InputParameters
LynxADHydroPoroMech::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Poromechanics coupling kernel.");
  return params;
}

LynxADHydroPoroMech::LynxADHydroPoroMech(const InputParameters & parameters)
  : ADKernel(parameters), _poro_mech(getADMaterialProperty<Real>("poro_mechanical"))
{
}

ADReal
LynxADHydroPoroMech::computeQpResidual()
{
  return _poro_mech[_qp] * _test[_i][_qp];
}