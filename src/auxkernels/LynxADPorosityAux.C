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

#include "LynxADPorosityAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LynxApp", LynxADPorosityAux);

InputParameters
LynxADPorosityAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Updating the porosity auxiliary variable.");
  params.addCoupledVar("fluid_pressure", 0.0, "The fluid pressure variable.");
  return params;
}

LynxADPorosityAux::LynxADPorosityAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _pf_dot(coupledDot("fluid_pressure")),
    _biot(getMaterialProperty<Real>("biot_coefficient")),
    _C_d(getMaterialProperty<Real>("bulk_compressibility")),
    _strain_increment(getADMaterialProperty<RankTwoTensor>("strain_increment")),
    // TO CORRECT with hasADMaterialProperty
    _has_viscous(hasMaterialProperty<RankTwoTensor>("viscous_strain_increment")),
    _viscous_strain_incr(
        _has_viscous ? &getADMaterialProperty<RankTwoTensor>("viscous_strain_increment") : nullptr),
    // TO CORRECT with hasADMaterialProperty
    _has_plastic(hasMaterialProperty<RankTwoTensor>("plastic_strain_increment")),
    _plastic_strain_incr(
        _has_viscous ? &getADMaterialProperty<RankTwoTensor>("plastic_strain_increment") : nullptr)
{
}

Real
LynxADPorosityAux::computeValue()
{
  Real ev_dot = computeEvDot();
  Real ev_in_dot = computeEvInDot();

  return _u_old[_qp] +
         (_biot[_qp] - _u[_qp]) *
             (_C_d[_qp] * (1.0 - _biot[_qp]) *
                  _pf_dot[_qp] +
              ev_dot) *
             _dt +
         (1.0 - _biot[_qp]) * ev_in_dot * _dt;
}

Real
LynxADPorosityAux::computeEvDot()
{
  ADRankTwoTensor e_dot = _strain_increment[_qp] / _dt;
  return MetaPhysicL::raw_value(e_dot.trace());
}

Real
LynxADPorosityAux::computeEvInDot()
{
  ADRankTwoTensor e_in_dot = ADRankTwoTensor();
  if (_has_viscous)
    e_in_dot += (*_viscous_strain_incr)[_qp] / _dt;
  if (_has_plastic)
    e_in_dot += (*_plastic_strain_incr)[_qp] / _dt;

  return MetaPhysicL::raw_value(e_in_dot.trace());
}
