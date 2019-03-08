/******************************************************************************/
/*                       LYNX, a MOOSE-based application                      */
/*                                                                            */
/*          Copyright (C) 2017 by Antoine B. Jacquey and Mauro Cacace         */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*    This program is free software: you can redistribute it and/or modify    */
/*    it under the terms of the GNU General Public License as published by    */
/*      the Free Software Foundation, either version 3 of the License, or     */
/*                     (at your option) any later version.                    */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program. If not, see <http://www.gnu.org/licenses/>     */
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
