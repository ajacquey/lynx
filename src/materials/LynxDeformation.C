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

#include "LynxDeformation.h"

registerMooseObject("LynxApp", LynxDeformation);

template <>
InputParameters
validParams<LynxDeformation>()
{
  InputParameters params = validParams<LynxDeformationBase>();
  params.addClassDescription("Class calculating the deformation of a material following a "
                             "viscoelasto-(visco)plastic rheology.");
  // Plastic parameters
  params.addParam<std::vector<Real>>(
      "friction_angle",
      "The friction angle of the material for the pressure-dependant part of the yield stress.");
  params.addParam<std::vector<Real>>("friction_angle_residual",
                                     "The residual friction angle of the material for the "
                                     "pressure-dependant part of the yield stress");
  params.addParam<std::vector<Real>>("cohesion",
                                     "The constant coefficient of the yield stress corresponding "
                                     "to the cohesion of the material.");
  params.addParam<std::vector<Real>>("cohesion_residual",
                                     "The residual of the constant coefficient of the yield stress "
                                     "corresponding to the cohesion of the material.");
  params.addParam<std::vector<Real>>(
      "dilation_angle",
      "The dilation angle of the material for the non-associative plastic update.");
  params.addParam<std::vector<Real>>(
      "internal_0", "The value of the plastic strain when hardening/softening begins.");
  params.addParam<std::vector<Real>>(
      "internal_limit", "The value of the plastic strain when hardening/softening ends.");
  params.addParam<std::vector<Real>>(
      "plastic_viscosity",
      "The reference viscosity in the generalized Duvaut-Lions viscoplastic formulation.");
  return params;
}

LynxDeformation::LynxDeformation(const InputParameters & parameters)
  : LynxDeformationBase(parameters),
    // Plastic parameters
    _friction_angle_0(isParamValid("friction_angle") ? getLynxParam<Real>("friction_angle")
                                                     : std::vector<Real>(_n_composition, 0.0)),
    _cohesion_0(isParamValid("cohesion") ? getLynxParam<Real>("cohesion")
                                         : std::vector<Real>(_n_composition, 0.0)),
    _friction_angle_res(isParamValid("friction_angle_residual")
                            ? getLynxParam<Real>("friction_angle_residual")
                            : _friction_angle_0),
    _cohesion_res(isParamValid("cohesion_residual") ? getLynxParam<Real>("cohesion_residual")
                                                    : _cohesion_0),
    _dilation_angle(isParamValid("dilation_angle") ? getLynxParam<Real>("dilation_angle")
                                                   : std::vector<Real>(_n_composition, 0.0)),
    _intnl_0(isParamValid("internal_0") ? getLynxParam<Real>("internal_0")
                                        : std::vector<Real>(_n_composition, 0.0)),
    _intnl_lim(isParamValid("internal_limit") ? getLynxParam<Real>("internal_limit")
                                              : std::vector<Real>(_n_composition, 1.0)),
    _one_on_plastic_eta(isParamValid("plastic_viscosity") ? getLynxParam<Real>("plastic_viscosity")
                                                          : std::vector<Real>(_n_composition, 0.0)),
    _has_hardening(((_friction_angle_0 != _friction_angle_res) || (_cohesion_0 != _cohesion_res)) &&
                   (_intnl_0 != _intnl_lim)),
    // Plastic properties
    _intnl(_has_hardening ? &declareProperty<Real>("plastic_intnl") : NULL),
    _intnl_old(_has_hardening ? &getMaterialPropertyOld<Real>("plastic_intnl") : NULL)
{
  _has_plasticity = (isParamValid("friction_angle") || isParamValid("cohesion"));

  // Grabing the inverse of plastic viscosity
  for (unsigned int i = 0; i < _n_composition; ++i)
    if (_one_on_plastic_eta[i] < 0.0)
      mooseError("LynxDeformation: 'plastic_viscosity' cannot be negative!");
    else if (_one_on_plastic_eta[i] > 0.0)
      _one_on_plastic_eta[i] = 1.0 / _one_on_plastic_eta[i];

  // Plastic structure
  if (_has_plasticity)
    _plasticity = new plasticity();
}

void
LynxDeformation::initQpStatefulProperties()
{
  if (_has_hardening)
    (*_intnl)[_qp] = 0.0;

  LynxDeformationBase::initQpStatefulProperties();
}

void
LynxDeformation::plasticCorrection(Real & pressure, RankTwoTensor & stress_dev)
{
  Real eqv_stress = std::sqrt(1.5) * stress_dev.L2norm();
  RankTwoTensor flow_dir = (eqv_stress != 0.0) ? stress_dev / eqv_stress : RankTwoTensor();

  computePlasticityProperties(pressure);

  _plastic_yield_function[_qp] = computePlasticityYield(pressure, eqv_stress);

  // Plastic correction
  if (_plastic_yield_function[_qp] > 0.0)
  {
    Real plastic_incr = plasticIncrement(pressure, eqv_stress);

    _plastic_strain_incr[_qp] = 1.5 * plastic_incr * flow_dir;
    _plastic_strain_incr[_qp].addIa(_plasticity->_beta * plastic_incr);
    stress_dev -= 3.0 * _G[_qp] * plastic_incr * flow_dir;
    pressure += _K[_qp] * _plasticity->_beta * plastic_incr;
    _elastic_strain[_qp] -= _plastic_strain_incr[_qp];
    // Update yield
    _plastic_yield_function[_qp] -= (3.0 * _G[_qp] + _plasticity->_H) * plastic_incr;
  }
}

Real
LynxDeformation::computePlasticityYield(const Real & pressure, const Real & eqv_stress)
{
  return eqv_stress - _plasticity->_alpha * pressure - _plasticity->_k;
}

Real
LynxDeformation::plasticIncrement(const Real & /*pressure*/, const Real & eqv_stress)
{
  Real eqv_p_strain_incr = 0.0;

  Real vp_correction = 1.0;
  eqv_p_strain_incr = _plastic_yield_function[_qp] / (3.0 * _G[_qp] + _plasticity->_H);

  // Viscoplastic correction
  if (_plasticity->_eta != 0.0)
  {
    vp_correction = (3.0 * _G[_qp] + _plasticity->_H) * _dt * _plasticity->_eta /
                    (1.0 + (3.0 * _G[_qp] + _plasticity->_H) * _dt * _plasticity->_eta);
    eqv_p_strain_incr *= vp_correction;
  }

  // Update internal parameter
  _plasticity->_intnl += eqv_p_strain_incr;
  if (_has_hardening)
    (*_intnl)[_qp] = _plasticity->_intnl;

  _stress_corr_p =
      (eqv_stress != 0.0) ? (eqv_stress - 3.0 * _G[_qp] * eqv_p_strain_incr) / eqv_stress : 1.0;
  _dp_dp_tr_p = 1.0 - _plasticity->_alpha * _plasticity->_beta * _K[_qp] * vp_correction /
                          (3.0 * _G[_qp] + _plasticity->_H);
  _dp_dq_tr_p = _plasticity->_beta * _K[_qp] * vp_correction / (3.0 * _G[_qp] + _plasticity->_H);
  _dq_dp_tr_p =
      _plasticity->_alpha * 3.0 * _G[_qp] * vp_correction / (3.0 * _G[_qp] + _plasticity->_H);
  _dq_dq_tr_p = 1.0 - 3.0 * _G[_qp] * vp_correction / (3.0 * _G[_qp] + _plasticity->_H);

  return eqv_p_strain_incr;
}

void
LynxDeformation::computePlasticityProperties(const Real & pressure)
{
  updatePlasticityParameters();
  if (_has_hardening)
    _plasticity->_intnl = (*_intnl_old)[_qp];

  if (_plasticity->_intnl < _plasticity->_intnl_0)
  {
    _plasticity->_alpha = _plasticity->_alpha_0;
    _plasticity->_k = _plasticity->_k_0;
    _plasticity->_H = 0.0;
    _plasticity->_H += _plasticity->_alpha * _plasticity->_beta * _K[_qp];
  }
  else if (_plasticity->_intnl >= _plasticity->_intnl_lim)
  {
    _plasticity->_alpha = _plasticity->_alpha_res;
    _plasticity->_k = _plasticity->_k_res;
    _plasticity->_H = 0.0;
    _plasticity->_H += _plasticity->_alpha * _plasticity->_beta * _K[_qp];
  }
  else
  {
    const Real x = (_plasticity->_intnl - _plasticity->_intnl_0) /
                   (_plasticity->_intnl_lim - _plasticity->_intnl_0);
    _plasticity->_alpha =
        _plasticity->_alpha_0 + (_plasticity->_alpha_res - _plasticity->_alpha_0) * x;
    _plasticity->_k = _plasticity->_k_0 + (_plasticity->_k_res - _plasticity->_k_0) * x;
    _plasticity->_H = (_plasticity->_alpha_res - _plasticity->_alpha_0) /
                          (_plasticity->_intnl_lim - _plasticity->_intnl_0) * pressure +
                      (_plasticity->_k_res - _plasticity->_k_0) /
                          (_plasticity->_intnl_lim - _plasticity->_intnl_0);
    _plasticity->_H +=
        _plasticity->_beta * _K[_qp] *
        (_plasticity->_alpha + (_plasticity->_alpha_res - _plasticity->_alpha_0) /
                                   (_plasticity->_intnl_lim - _plasticity->_intnl_0));
  }
}

void
LynxDeformation::updatePlasticityParameters()
{
  Real alpha_0 =
      std::sqrt(3.0) * std::sin(averageProperty(_friction_angle_0) * libMesh::pi / 180.0);
  Real alpha_res =
      std::sqrt(3.0) * std::sin(averageProperty(_friction_angle_res) * libMesh::pi / 180.0);
  Real k_0 = std::sqrt(3.0) * averageProperty(_cohesion_0) *
             std::cos(averageProperty(_friction_angle_0) * libMesh::pi / 180.0);
  Real k_res = std::sqrt(3.0) * averageProperty(_cohesion_res) *
               std::cos(averageProperty(_friction_angle_res) * libMesh::pi / 180.0);
  Real beta = std::sqrt(3.0) * std::sin(averageProperty(_dilation_angle) * libMesh::pi / 180.0);
  _plasticity->fill(alpha_0,
                    alpha_res,
                    k_0,
                    k_res,
                    beta,
                    averageProperty(_intnl_0),
                    averageProperty(_intnl_lim),
                    averageProperty(_one_on_plastic_eta));
}
