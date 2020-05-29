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

#include "LynxADPlasticModel.h"

registerMooseObject("LynxApp", LynxADPlasticModel);

InputParameters
LynxADPlasticModel::validParams()
{
  InputParameters params = LynxADMaterialBase::validParams();
  params.addClassDescription("Base class for the plastic correction of a elasto-plastic rheology.");
  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");
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

LynxADPlasticModel::LynxADPlasticModel(const InputParameters & parameters)
  : LynxADMaterialBase(parameters),
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
    _intnl_param_0(isParamValid("internal_0") ? getLynxParam<Real>("internal_0")
                                              : std::vector<Real>(_n_composition, 0.0)),
    _intnl_param_lim(isParamValid("internal_limit") ? getLynxParam<Real>("internal_limit")
                                                    : std::vector<Real>(_n_composition, 1.0)),
    _one_on_plastic_eta(isParamValid("plastic_viscosity") ? getLynxParam<Real>("plastic_viscosity")
                                                          : std::vector<Real>(_n_composition, 0.0)),
    _has_hardening(((_friction_angle_0 != _friction_angle_res) || (_cohesion_0 != _cohesion_res)) &&
                   (_intnl_param_0 != _intnl_param_lim)),
    // PLastic properties
    _plastic_yield_function(declareADProperty<Real>("plastic_yield_function")),
    _plastic_strain_incr(declareADProperty<RankTwoTensor>("plastic_strain_increment")),
    _intnl(_has_hardening ? &declareADProperty<Real>("plastic_intnl") : nullptr),
    _intnl_old(_has_hardening ? &getMaterialPropertyOld<Real>("plastic_intnl") : nullptr)
{
  // Grabing the inverse of plastic viscosity
  for (unsigned int i = 0; i < _n_composition; ++i)
    if (_one_on_plastic_eta[i] < 0.0)
      mooseError("LynxPlasticModel: 'plastic_viscosity' cannot be negative!");
    else if (_one_on_plastic_eta[i] > 0.0)
      _one_on_plastic_eta[i] = 1.0 / _one_on_plastic_eta[i];
}

void
LynxADPlasticModel::initQpStatefulProperties()
{
  if (_has_hardening)
    (*_intnl)[_qp] = 0.0;
}

void
LynxADPlasticModel::setQp(unsigned int qp)
{
  _qp = qp;
}

void
LynxADPlasticModel::plasticUpdate(ADRankTwoTensor & stress_dev,
                                  ADReal & pressure,
                                  const Real & G,
                                  const Real & K,
                                  ADRankTwoTensor & elastic_strain_incr)
{
  const ADReal eqv_stress = std::sqrt(1.5) * stress_dev.L2norm();
  const ADRankTwoTensor flow_dir =
      (eqv_stress != 0.0) ? stress_dev / eqv_stress : ADRankTwoTensor();

  ADReal delta_e_eqv = plasticIncrement(eqv_stress, pressure, G, K);

  _plastic_strain_incr[_qp] = 1.5 * delta_e_eqv * flow_dir;
  _plastic_strain_incr[_qp].addIa(_beta * delta_e_eqv / 3.0);
  stress_dev -= 3.0 * G * delta_e_eqv * flow_dir;
  pressure += K * _beta * delta_e_eqv;
  elastic_strain_incr -= _plastic_strain_incr[_qp];
  // Update yield function
  initPlasticParameters(pressure, K);
  _plastic_yield_function[_qp] =
      plasticYieldFunction(std::sqrt(1.5) * stress_dev.L2norm(), pressure);
}

ADReal
LynxADPlasticModel::plasticIncrement(const ADReal & eqv_stress,
                                     const ADReal & pressure,
                                     const Real G,
                                     const Real K)
{
  // Initialize hardening
  if (_has_hardening)
    (*_intnl)[_qp] = (*_intnl_old)[_qp];

  // Map plastic parameters
  initPlasticParameters(pressure, K);

  // If no friction and cohesion are provided, there is no plastic update
  if ((_alpha == 0.0 && _k == 0.0) || G == 0.0)
  {
    _plastic_yield_function[_qp] = -1.0;
    return 0.0;
  }

  // Yield function
  _plastic_yield_function[_qp] = plasticYieldFunction(eqv_stress, pressure);

  ADReal delta_e_eqv = 0.0;
  if (_plastic_yield_function[_qp] <= 0.0) // Elastic
    return 0.0;

  // Plastic step
  delta_e_eqv = _plastic_yield_function[_qp] / (3.0 * G + _H);

  // Visco-plastic correction
  // Viscoplastic correction
  if (_one_on_eta != 0.0)
  {
    ADReal vp_correction =
        (3.0 * G + _H) * _dt * _one_on_eta / (1.0 + (3.0 * G + _H) * _dt * _one_on_eta);
    delta_e_eqv *= vp_correction;
  }

  // Update internal parameter
  if (_has_hardening)
    (*_intnl)[_qp] += delta_e_eqv;

  return delta_e_eqv;
}

void
LynxADPlasticModel::initPlasticParameters(const ADReal & pressure, const Real & K)
{
  _alpha_0 = std::sqrt(3.0) * std::sin(averageProperty(_friction_angle_0) * libMesh::pi / 180.0);
  _alpha_res =
      std::sqrt(3.0) * std::sin(averageProperty(_friction_angle_res) * libMesh::pi / 180.0);
  _k_0 = std::sqrt(3.0) * averageProperty(_cohesion_0) *
         std::cos(averageProperty(_friction_angle_0) * libMesh::pi / 180.0);
  _k_res = std::sqrt(3.0) * averageProperty(_cohesion_res) *
           std::cos(averageProperty(_friction_angle_res) * libMesh::pi / 180.0);
  _beta = std::sqrt(3.0) * std::sin(averageProperty(_dilation_angle) * libMesh::pi / 180.0);
  _intnl_0 = averageProperty(_intnl_param_0);
  _intnl_lim = averageProperty(_intnl_param_lim);
  _one_on_eta = averageProperty(_one_on_plastic_eta);

  _alpha = _alpha_0;
  _k = _k_0;
  _H = _alpha * _beta * K;

  if (_has_hardening)
  {
    if ((*_intnl)[_qp] >= _intnl_lim)
    {
      _alpha = _alpha_res;
      _k = _k_res;
      _H = _alpha * _beta * K;
    }
    else
    {
      const ADReal x = ((*_intnl)[_qp] - _intnl_0) / (_intnl_lim - _intnl_0);
      _alpha = _alpha_0 + (_alpha_res - _alpha_0) * x;
      _k = _k_0 + (_k_res - _k_0) * x;
      _H = (_alpha_res - _alpha_0) / (_intnl_lim - _intnl_0) * pressure +
           (_k_res - _k_0) / (_intnl_lim - _intnl_0);
      _H += _beta * K * (_alpha + (_alpha_res - _alpha_0) / (_intnl_lim - _intnl_0));
    }
  }
}

ADReal
LynxADPlasticModel::plasticYieldFunction(const ADReal & eqv_stress, const ADReal & pressure)
{
  return eqv_stress - _alpha * pressure - _k;
}