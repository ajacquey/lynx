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

#include "LynxHoldStressBC.h"
#include "Function.h"

registerMooseObject("LynxApp", LynxHoldStressBC);

InputParameters
LynxHoldStressBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription("Hold the stress on a given boundary in a given direction.");
  params.addCoupledVar("fluid_pressure", "The fluid pressure variable.");
  params.addRequiredParam<unsigned int>("component", "The component for the pressure.");
  // Elastic moduli parameters
  params.addRangeCheckedParam<Real>(
      "bulk_modulus", "bulk_modulus >= 0.0", "The drained bulk modulus of the material.");
  params.addRangeCheckedParam<Real>(
      "shear_modulus", "shear_modulus >= 0.0", "The shear modulus of the material.");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

LynxHoldStressBC::LynxHoldStressBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _coupled_pf(isCoupled("fluid_pressure")),
    _pf(_coupled_pf ? coupledValueOld("fluid_pressure") : _zero),
    _component(getParam<unsigned int>("component")),
    _bulk_modulus(getParam<Real>("bulk_modulus")),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _elastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("elastic_strain")),
    _biot_coeff(getMaterialProperty<Real>("biot_coefficient"))
{
  if (_component > 2)
    mooseError("Invalid component given for ", name(), ": ", _component, ".\n");
}

Real
LynxHoldStressBC::computeQpResidual()
{
  Real vol_strain = _elastic_strain_old[_qp].trace();
  RankTwoTensor dev_strain = _elastic_strain_old[_qp].deviatoric();

  RankTwoTensor stress = 2.0 * _shear_modulus * dev_strain;
  stress.addIa(_bulk_modulus * vol_strain);

  Real value = -(stress(_component, _component) - _biot_coeff[_qp] * _pf[_qp]);

  return value * (_normals[_qp](_component) * _test[_i][_qp]);
}
