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

#include "LynxWinklerBC.h"
#include "Function.h"

registerMooseObject("LynxApp", LynxWinklerBC);

template <>
InputParameters
validParams<LynxWinklerBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addClassDescription(
      "Applies a Winkler foundation BC on a given boundary in a given direction.");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system.");
  params.addRequiredParam<unsigned int>(
      "component", "The component of the displacement orthogonal to the boundary.");
  params.addParam<Real>("value", 0.0, "Value of the reference pressure.");
  params.addParam<FunctionName>("function", "The function that describes the reference pressure.");
  params.addParam<Real>("external_density", 0.0, "The density of the external material.");
  params.addParam<Real>("gravity_acceleration", 9.81, "The magnitude of the gravity acceleration.");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

LynxWinklerBC::LynxWinklerBC(const InputParameters & parameters)
  : DerivativeMaterialInterface<IntegratedBC>(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp(3),
    _component(getParam<unsigned int>("component")),
    _value(getParam<Real>("value")),
    _function(isParamValid("function") ? &getFunction("function") : NULL),
    _rho_ext(getParam<Real>("external_density")),
    _g(getParam<Real>("gravity_acceleration")),
    _rho_b(getMaterialProperty<Real>("bulk_density"))
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError(
        "The number of variables supplied in 'displacements' must match the mesh dimension.");

  // Fetch coupled variables
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp[i] = &coupledValue("displacements", i);
  // Set unused dimensions to zero
  for (unsigned i = _ndisp; i < 3; ++i)
    _disp[i] = &_zero;

  if (_component > 2)
    mooseError("Invalid component given for ", name(), ": ", _component, ".\n");
}

Real
LynxWinklerBC::computeQpResidual()
{
  Real value = _value;

  if (_function)
    value = _function->value(_t, _q_point[_qp]);

  // Correct pressure for external density
  if (_rho_ext != 0.0)
    value += (_rho_ext - _rho_b[_qp]) * _g * (*_disp[_component])[_qp];
  return value * (_normals[_qp](_component) * _test[_i][_qp]);
}
