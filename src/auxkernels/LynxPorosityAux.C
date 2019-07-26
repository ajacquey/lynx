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

#include "LynxPorosityAux.h"

registerMooseObject("LynxApp", LynxPorosityAux);

template <>
InputParameters
validParams<LynxPorosityAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Updating the porosity auxiliary variable.");
  params.addCoupledVar("fluid_pressure", "The fluid pressure variable.");
  return params;
}

LynxPorosityAux::LynxPorosityAux(const InputParameters & parameters)
  : DerivativeMaterialInterface<AuxKernel>(parameters),
    _coupled_pf(isCoupled("fluid_pressure")),
    _pf_dot(_coupled_pf ? coupledDot("fluid_pressure") : _zero),
    _biot(getDefaultMaterialProperty<Real>("biot_coefficient")),
    _C_d(getDefaultMaterialProperty<Real>("bulk_compressibility")),
    _strain_increment(getDefaultMaterialProperty<RankTwoTensor>("strain_increment")),
    _inelastic_strain(getDefaultMaterialProperty<RankTwoTensor>("inelastic_strain")),
    _inelastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("inelastic_strain"))
{
}

Real
LynxPorosityAux::computeValue()
{
  Real ev_dot = computeEvDot();
  Real ev_in_dot = computeEvInDot();

  return _u_old[_qp] +
         (_biot[_qp] - _u[_qp]) * (_C_d[_qp] * (1.0 - _biot[_qp]) * _pf_dot[_qp] + ev_dot) * _dt +
         (1.0 - _biot[_qp]) * ev_in_dot * _dt;
}

Real
LynxPorosityAux::computeEvDot()
{
  RankTwoTensor e_dot = _strain_increment[_qp] / _dt;
  return e_dot.trace();
}

Real
LynxPorosityAux::computeEvInDot()
{
  RankTwoTensor e_in_dot = (_inelastic_strain[_qp] - _inelastic_strain_old[_qp]) / _dt;
  return e_in_dot.trace();
}
