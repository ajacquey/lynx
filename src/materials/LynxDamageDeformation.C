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

#include "LynxDamageDeformation.h"
#include "libmesh/utility.h"

registerMooseObject("LynxApp", LynxDamageDeformation);

template <>
InputParameters
validParams<LynxDamageDeformation>()
{
  InputParameters params = validParams<LynxDeformationBase>();
  params.addClassDescription("Class calculating the deformation of a material following a "
                             "damage viscoelasto-(visco)plastic rheology.");
  // Coupled variables
  params.addCoupledVar("damage", "The damage variable.");
  params.addCoupledVar("porosity", "The porosity variable.");
  // Elastic moduli parameters
  params.addParam<std::vector<Real>>("damage_modulus",
                                     "The third elastic damage modulus for the damage rheology.");
  // Damage-Plasticity parameters
  params.addParam<std::vector<Real>>(
      "friction_angle",
      "The friction angle of the material for the pressure-dependant part of the yield stress.");
  params.addParam<std::vector<Real>>(
      "porous_coeff",
      "The porous coefficient controlling the capped yield in the damage rheology.");
  params.addParam<std::vector<Real>>(
      "porous_coeff_linear", "The linear dependency of the porous coefficient to porosity.");
  params.addParam<std::vector<Real>>(
      "plastic_viscosity",
      "The reference viscosity in the generalized Duvaut-Lions viscoplastic formulation.");
  params.addParam<std::vector<Real>>("damage_viscosity",
                                     "The reference viscosity in the damage rheology.");
  return params;
}

LynxDamageDeformation::LynxDamageDeformation(const InputParameters & parameters)
  : LynxDeformationBase(parameters),
    // Coupled variables
    _coupled_dam(isCoupled("damage")),
    _damage(_coupled_dam ? coupledValue("damage") : _zero),
    _damage_old(_coupled_dam ? coupledValueOld("damage") : _zero),
    _coupled_phi(isCoupled("porosity")),
    _porosity(_coupled_phi ? coupledValue("porosity") : _zero),
    // Elastic moduli parameters
    _damage_modulus(isParamValid("damage_modulus") ? getLynxParam<Real>("damage_modulus")
                                                   : std::vector<Real>(_n_composition, 0.0)),

    // Damage-Plasticity parameters
    _friction_angle(isParamValid("friction_angle") ? getLynxParam<Real>("friction_angle")
                                                   : std::vector<Real>(_n_composition, 0.0)),
    _porous_coeff(isParamValid("porous_coeff") ? getLynxParam<Real>("porous_coeff")
                                               : std::vector<Real>(_n_composition, 0.0)),
    _porous_coeff_linear(isParamValid("porous_coeff_linear")
                             ? getLynxParam<Real>("porous_coeff_linear")
                             : std::vector<Real>(_n_composition, 0.0)),
    _one_on_plastic_eta(isParamValid("plastic_viscosity") ? getLynxParam<Real>("plastic_viscosity")
                                                          : std::vector<Real>(_n_composition, 0.0)),
    _one_on_damage_eta(isParamValid("damage_viscosity") ? getLynxParam<Real>("damage_viscosity")
                                                        : std::vector<Real>(_n_composition, 0.0)),
    // Strain properties
    _elastic_strain(declareProperty<RankTwoTensor>("elastic_strain")),
    _elastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("elastic_strain"))
{
  _has_plasticity = isParamValid("friction_angle");

  // Grabing the inverse of plastic and damage viscosities
  for (unsigned int i = 0; i < _n_composition; ++i)
  {
    if (_one_on_plastic_eta[i] < 0.0)
      mooseError("LynxDamageDeformation: 'plastic_viscosity' cannot be negative!");
    else if (_one_on_plastic_eta[i] > 0.0)
      _one_on_plastic_eta[i] = 1.0 / _one_on_plastic_eta[i];

    if (_one_on_damage_eta[i] < 0.0)
      mooseError("LynxDamageDeformation: 'damage_viscosity' cannot be negative!");
    else if (_one_on_damage_eta[i] > 0.0)
      _one_on_damage_eta[i] = 1.0 / _one_on_damage_eta[i];
  }

  // Damage-Plasticity structure
  if (_has_plasticity)
    _damage_plasticity = new damage_plasticity();
}

void
LynxDamageDeformation::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
  LynxDeformationBase::initQpStatefulProperties();
}

void
LynxDamageDeformation::initializeQpDeformation()
{
  // Need to get _xi0
  initializeDamageParameters();

  LynxDeformationBase::initializeQpDeformation();
}

void
LynxDamageDeformation::elasticModuli()
{
  // Strain ratio
  const Real xi = strainRatio(_elastic_strain_old[_qp]);

  // Bulk modulus
  _K[_qp] = _has_bulk_modulus ? averageProperty(_bulk_modulus) -
                                    2.0 / 3.0 * _damage_old[_qp] * _damage_plasticity->_gamma *
                                        (xi / 2.0 - _damage_plasticity->_xi0)
                              : 0.0;

  // Shear modulus
  _G[_qp] = _has_shear_modulus
                ? averageProperty(_shear_modulus) - _damage_old[_qp] * _damage_plasticity->_gamma *
                                                        (xi / 2.0 - _damage_plasticity->_xi0)
                : 0.0;
}

void
LynxDamageDeformation::plasticCorrection(Real & pressure, RankTwoTensor & stress_dev)
{
  Real eqv_stress = std::sqrt(1.5) * stress_dev.L2norm();
  RankTwoTensor flow_dir = (eqv_stress != 0.0) ? stress_dev / eqv_stress : RankTwoTensor();

  computeDamageProperties(pressure, eqv_stress);

  // Check wether the yield is capped
  if (_damage_plasticity->_p_cr == 0.0)
  {
    // No capped yield - only deviatoric correction
    _plastic_yield_function[_qp] = computePlasticityYield(pressure, eqv_stress);

    // Plastic correction
    if (_plastic_yield_function[_qp] > 0.0)
    {
      Real plastic_incr = plasticIncrement(pressure, eqv_stress);

      _plastic_strain_incr[_qp] = 1.5 * plastic_incr * flow_dir;
      stress_dev -= 3.0 * _G[_qp] * plastic_incr * flow_dir;
      // Update yield
      _plastic_yield_function[_qp] -= 3.0 * _G[_qp] * plastic_incr;
    }
  }
  else
  {
    // Capped yield
    Real yield_squared = computeConvexPlasticityYield2(pressure, eqv_stress);
    _plastic_yield_function[_qp] =
        (yield_squared > 0.0) ? std::sqrt(yield_squared) : -std::sqrt(-yield_squared);

    if (_plastic_yield_function[_qp] > 0.0)
    {
      Real vol_plastic_incr = 0.0, eqv_plastic_incr = 0.0;
      Real rho_0 = convexPlasticIncrement(vol_plastic_incr, eqv_plastic_incr);

      _plastic_strain_incr[_qp] = 1.5 * eqv_plastic_incr * flow_dir;
      _plastic_strain_incr[_qp].addIa(-vol_plastic_incr / 3.0);

      pressure -= _K[_qp] * vol_plastic_incr;
      stress_dev -= 3.0 * _G[_qp] * eqv_plastic_incr * flow_dir;
      eqv_stress = std::sqrt(1.5) * stress_dev.L2norm();
      // Update yield
      _plastic_yield_function[_qp] =
          std::sqrt(Utility::pow<2>(pressure - _damage_plasticity->_p_r) +
                    Utility::pow<2>(eqv_stress - _damage_plasticity->_q_r)) -
          rho_0;
    }
  }

  // Update elastic strain
  _elastic_strain[_qp] = spinRotation(_elastic_strain_old[_qp].deviatoric());
  _elastic_strain[_qp].addIa(_elastic_strain_old[_qp].trace() / 3.0);

  _elastic_strain[_qp] +=
      _strain_increment[_qp] - _viscous_strain_incr[_qp] - _plastic_strain_incr[_qp];
}

void
LynxDamageDeformation::damageCorrection()
{
}

Real
LynxDamageDeformation::computePlasticityYield(const Real & pressure, const Real & eqv_stress)
{
  return eqv_stress - _damage_plasticity->_alpha0 * pressure;
}

Real
LynxDamageDeformation::plasticIncrement(const Real & /*pressure*/, const Real & eqv_stress)
{
  Real eqv_p_strain_incr = 0.0;

  Real vp_correction = 1.0;
  eqv_p_strain_incr = _plastic_yield_function[_qp] / (3.0 * _G[_qp]);

  // Viscoplastic correction
  if (_damage_plasticity->_eta_p != 0.0)
  {
    vp_correction = (3.0 * _G[_qp]) * _dt * _damage_plasticity->_eta_p /
                    (1.0 + (3.0 * _G[_qp]) * _dt * _damage_plasticity->_eta_p);
    eqv_p_strain_incr *= vp_correction;
  }

  // Update yield
  _plastic_yield_function[_qp] -= (3.0 * _G[_qp]) * eqv_p_strain_incr;

  _stress_corr_p =
      (eqv_stress != 0.0) ? (eqv_stress - 3.0 * _G[_qp] * eqv_p_strain_incr) / eqv_stress : 1.0;
  _dp_dp_tr_p = 1.0;
  _dp_dq_tr_p = 0.0;
  _dq_dp_tr_p = _damage_plasticity->_alpha0 * vp_correction;
  _dq_dq_tr_p = 1.0 - vp_correction;

  return eqv_p_strain_incr;
}

Real
LynxDamageDeformation::computeConvexPlasticityYield2(const Real & pressure, const Real & eqv_stress)
{
  if ((pressure == 0.0) && (eqv_stress == 0.0))
    return 0.0;

  updateDamageConvexParameters(pressure, eqv_stress);

  return Utility::pow<2>(eqv_stress) - _damage_plasticity->_alpha2 * Utility::pow<2>(pressure) *
                                           ((0.0 < pressure) - (pressure < 0.0));
}

Real
LynxDamageDeformation::computeConvexPlasticityYield2(const Real & rho)
{
  Real p = _damage_plasticity->_p_r + rho * _damage_plasticity->_rp;
  Real q = _damage_plasticity->_q_r + rho * _damage_plasticity->_rq;

  return computeConvexPlasticityYield2(p, q);
}

Real
LynxDamageDeformation::convexPlasticIncrement(Real & vol_plastic_incr, Real & eqv_plastic_incr)
{
  Real vp_correction = 1.0;

  // Find projection on convex yield
  Real rho_neg = 0.0;
  Real rho_pos = _damage_plasticity->_rho_tr;
  Real rho_0 = getConvexProjection(rho_neg, rho_pos);

  // Update projections
  Real p_0 = _damage_plasticity->_p_r + rho_0 * _damage_plasticity->_rp;
  Real q_0 = _damage_plasticity->_q_r + rho_0 * _damage_plasticity->_rq;

  // Plastic strain
  vol_plastic_incr = (_damage_plasticity->_p_tr - p_0) / _K[_qp];
  eqv_plastic_incr = (_damage_plasticity->_q_tr - q_0) / (3.0 * _G[_qp]);

  // Viscoplastic correction
  if (_damage_plasticity->_eta_p != 0.0)
  {
    Real F_tr = _damage_plasticity->_rho_tr - rho_0;
    Real delta_p = std::sqrt(Utility::pow<2>(vol_plastic_incr) + Utility::pow<2>(eqv_plastic_incr));

    vp_correction = (F_tr * _damage_plasticity->_eta_p * _dt / delta_p) /
                    (1.0 + F_tr * _damage_plasticity->_eta_p * _dt / delta_p);

    vol_plastic_incr *= vp_correction;
    eqv_plastic_incr *= vp_correction;
  }

  // Tangent operator
  _stress_corr_p =
      (_damage_plasticity->_q_tr != 0.0)
          ? _damage_plasticity->_q_tr - 3.0 * _G[_qp] * eqv_plastic_incr / _damage_plasticity->_q_tr
          : 1.0;
  Real dF2_dp0 = dConvexPlasticYield2_dp(p_0, q_0) / dConvexPlasticYield2(rho_0);
  Real dF2_dq0 = dConvexPlasticYield2_dq(p_0, q_0) / dConvexPlasticYield2(rho_0);

  // _dp_dp_tr_p =
  //     _dp_r_dp_tr + (1.0 - vp_correction) + vp_correction * rho_0 / _rho_tr * dF2_dq0 * _rq;
  // _dp_dq_tr_p = _dp_r_dq_tr - vp_correction * rho_0 / _rho_tr * dF2_dq0 * _rp;
  // _dq_dp_tr_p = -vp_correction * rho_0 / _rho_tr * dF2_dp0 * _rq;
  // _dq_dq_tr_p = (1.0 - vp_correction) + vp_correction * rho_0 / _rho_tr * dF2_dp0 * _rp;
  _dp_dp_tr_p = (1.0 - vp_correction) + vp_correction * rho_0 / _damage_plasticity->_rho_tr *
                                            dF2_dq0 * _damage_plasticity->_rq;
  _dp_dq_tr_p =
      -vp_correction * rho_0 / _damage_plasticity->_rho_tr * dF2_dq0 * _damage_plasticity->_rp;
  _dq_dp_tr_p =
      -vp_correction * rho_0 / _damage_plasticity->_rho_tr * dF2_dp0 * _damage_plasticity->_rq;
  _dq_dq_tr_p = (1.0 - vp_correction) + vp_correction * rho_0 / _damage_plasticity->_rho_tr *
                                            dF2_dp0 * _damage_plasticity->_rp;

  return rho_0;
}

void
LynxDamageDeformation::computeDamageProperties(const Real & pressure, const Real & eqv_stress)
{
  updateDamageParameters();

  _damage_plasticity->_alpha0 = std::sqrt(2.0) * _G[_qp] / _K[_qp] *
                                std::sqrt((3.0 - Utility::pow<2>(_damage_plasticity->_xi0)) /
                                          Utility::pow<2>(_damage_plasticity->_xi0));

  // Properties for convex yield
  if (_damage_plasticity->_p_cr != 0.0)
  {
    // Trial pressure - deviatoric stress
    _damage_plasticity->_p_tr = pressure;
    _damage_plasticity->_q_tr = eqv_stress;

    // Reference point for convex yield
    _damage_plasticity->_p_r = convexReferencePressure();
    _damage_plasticity->_q_r = 0.0;

    // Trial distance
    _damage_plasticity->_rho_tr = std::sqrt(Utility::pow<2>(pressure - _damage_plasticity->_p_r) +
                                            Utility::pow<2>(eqv_stress - _damage_plasticity->_q_r));

    // Projection direction
    _damage_plasticity->_rp =
        (_damage_plasticity->_p_tr - _damage_plasticity->_p_r) / _damage_plasticity->_rho_tr;
    _damage_plasticity->_rq =
        (_damage_plasticity->_q_tr - _damage_plasticity->_q_r) / _damage_plasticity->_rho_tr;
  }
}

void
LynxDamageDeformation::updateDamageParameters()
{
  Real sin_phi = std::sin(averageProperty(_friction_angle) * libMesh::pi / 180.0);
  Real xi0 = -std::sqrt(3.0) / std::sqrt(1.0 + 1.5 * Utility::pow<2>(_K[_qp] / _G[_qp] * sin_phi));
  Real porous_coeff =
      averageProperty(_porous_coeff) + averageProperty(_porous_coeff_linear) * _porosity[_qp];
  Real p_cr = _K[_qp] * (std::sqrt(3.0) + xi0) / (3.0 * porous_coeff);
  _damage_plasticity->fill(xi0,
                           averageProperty(_damage_modulus),
                           p_cr,
                           averageProperty(_one_on_plastic_eta),
                           averageProperty(_one_on_damage_eta));
}

void
LynxDamageDeformation::initializeDamageParameters()
{
  Real sin_phi = std::sin(averageProperty(_friction_angle) * libMesh::pi / 180.0);
  Real xi0 = -std::sqrt(3.0) /
             std::sqrt(1.0 + 1.5 * Utility::pow<2>(averageProperty(_bulk_modulus) /
                                                   averageProperty(_shear_modulus) * sin_phi));
  Real porous_coeff =
      averageProperty(_porous_coeff) + averageProperty(_porous_coeff_linear) * _porosity[_qp];
  Real p_cr = _K[_qp] * (std::sqrt(3.0) + xi0) / (3.0 * porous_coeff);
  _damage_plasticity->fill(xi0,
                           averageProperty(_damage_modulus),
                           p_cr,
                           averageProperty(_one_on_plastic_eta),
                           averageProperty(_one_on_damage_eta));
}

void
LynxDamageDeformation::updateDamageConvexParameters(const Real & pressure, const Real & eqv_stress)
{
  Real pq_term =
      (pressure != 0.0)
          ? 1.0 / (1.0 + 0.5 * Utility::pow<2>(_K[_qp] * eqv_stress / (_G[_qp] * pressure)))
          : 0.0;

  _damage_plasticity->_xi_cr =
      (pressure > 0.0)
          ? _damage_plasticity->_xi0 - (std::sqrt(3.0) + _damage_plasticity->_xi0) * pressure /
                                           _damage_plasticity->_p_cr * pq_term
          : _damage_plasticity->_xi0;

  _damage_plasticity->_alpha2 = (_damage_plasticity->_xi_cr != 0.0)
                                    ? 2.0 * Utility::pow<2>(_G[_qp] / _K[_qp]) *
                                          (3.0 - Utility::pow<2>(_damage_plasticity->_xi_cr)) /
                                          Utility::pow<2>(_damage_plasticity->_xi_cr)
                                    : 0.0;
  _damage_plasticity->_dxi_cr_dp = (pressure != 0.0)
                                       ? -(_damage_plasticity->_xi0 - _damage_plasticity->_xi_cr) /
                                             pressure * (3.0 - 2.0 * pq_term)
                                       : 0.0;

  _damage_plasticity->_dxi_cr_dq =
      (eqv_stress != 0.0) ? 2.0 * (_damage_plasticity->_xi0 - _damage_plasticity->_xi_cr) /
                                eqv_stress * (1.0 - pq_term)
                          : 0.0;

  _damage_plasticity->_dmu2_dxi_cr =
      (_damage_plasticity->_xi_cr != 0.0)
          ? -12.0 * Utility::pow<2>(_G[_qp] / _K[_qp]) / Utility::pow<3>(_damage_plasticity->_xi_cr)
          : 0.0;
}

Real
LynxDamageDeformation::convexReferencePressure()
{
  Real p_proj =
      (1.0 - Utility::pow<2>(_damage_plasticity->_alpha0) * _K[_qp] /
                 (3.0 * _G[_qp] + Utility::pow<2>(_damage_plasticity->_alpha0) * _K[_qp])) *
          _damage_plasticity->_p_tr +
      _damage_plasticity->_alpha0 * _K[_qp] * _damage_plasticity->_q_tr /
          (3.0 * _G[_qp] + Utility::pow<2>(_damage_plasticity->_alpha0) * _K[_qp]);

  return ((p_proj > 0.0) && (p_proj < _damage_plasticity->_p_cr))
             ? p_proj
             : 2.0 * _damage_plasticity->_p_cr / 3.0;
}

Real
LynxDamageDeformation::dConvexPlasticYield2(const Real & rho)
{
  Real p = _damage_plasticity->_p_r + rho * _damage_plasticity->_rp;
  Real q = _damage_plasticity->_q_r + rho * _damage_plasticity->_rq;

  return dConvexPlasticYield2_dp(p, q) * _damage_plasticity->_rp +
         dConvexPlasticYield2_dq(p, q) * _damage_plasticity->_rq;
}

Real
LynxDamageDeformation::dConvexPlasticYield2_dp(const Real & pressure, const Real & eqv_stress)
{
  if ((pressure == 0.0) && (eqv_stress == 0.0))
    return 0.0;

  updateDamageConvexParameters(pressure, eqv_stress);

  return -(_damage_plasticity->_dmu2_dxi_cr * _damage_plasticity->_dxi_cr_dp * pressure +
           2.0 * _damage_plasticity->_alpha2) *
         pressure * ((0.0 < pressure) - (pressure < 0.0));
}

Real
LynxDamageDeformation::dConvexPlasticYield2_dq(const Real & pressure, const Real & eqv_stress)
{
  if ((pressure == 0.0) && (eqv_stress == 0.0))
    return 0.0;

  updateDamageConvexParameters(pressure, eqv_stress);

  return 2.0 * eqv_stress -
         _damage_plasticity->_dmu2_dxi_cr * _damage_plasticity->_dxi_cr_dq *
             Utility::pow<2>(pressure); // *
                                        // ((0.0 < pressure - _p_k) - (pressure - _p_k < 0.0));
}

Real
LynxDamageDeformation::getConvexProjection(const Real & x1, const Real & x2)
{
  // Machine floating-point precision
  const Real eps = std::numeric_limits<Real>::epsilon();
  // Declare some stuff here
  Real a = x1, b = x2, c = x2, d, e, fa = computeConvexPlasticityYield2(a),
       fb = computeConvexPlasticityYield2(b), fc, p, q, r, s, tol1, xm;
  // Chek if x1 and x2 bracket a root
  if (fa * fb > 0.0) // fa and fb are of the same sign
    throw MooseException("LynxDamageDeformation: in Brent's method, the points x1 and x2 must "
                         "bracket a root of the function!\n");

  fc = fb;
  for (unsigned int iter = 0; iter < _itmax; ++iter)
  {
    if (fb * fc > 0.0) // fb and fc are of the same sign
    {
      // Rename a, b, c and adjust bounding interval d
      c = a;
      fc = fa;
      e = d = b - a;
    }
    if (std::abs(fc) < std::abs(fb))
    {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    // Convergence check
    tol1 = 2.0 * eps * std::abs(b) + 0.5 * _tol;
    xm = 0.5 * (c - b);
    if (std::abs(xm) <= tol1 || fb == 0.0) // Exit of the algorithm
      return b;
    if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb))
    {
      s = fb / fa;
      if (a == c)
      {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      }
      else
      {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      // Check wether in bounds
      if (p > 0.0)
        q = -q;
      p = std::abs(p);
      Real min1 = 3.0 * xm * q - std::abs(tol1 * q);
      Real min2 = std::abs(e * q);
      if (2.0 * p < (min1 < min2 ? min1 : min2))
      {
        // Accept interpolation
        e = d;
        d = p / q;
      }
      else
      {
        // Interpolation failed, use bisection
        d = xm;
        e = d;
      }
    }
    else
    {
      // Bounds decreasing too slowly, use bisection
      d = xm;
      e = d;
    }
    // Move last guess to a
    a = b;
    fa = fb;
    // Evaluate new trial root
    if (std::abs(d) > tol1)
      b += d;
    else
      b += (xm < 0 ? -std::abs(tol1) : std::abs(tol1));
    fb = computeConvexPlasticityYield2(b);
  }
  throw MooseException(
      "LynxDamageDeformation: maximum number of iterations exceeded in solveBrent!");
}

Real
LynxDamageDeformation::strainRatio(const RankTwoTensor & elastic_strain)
{
  const Real strain_v = elastic_strain.trace();
  const Real strain_norm = elastic_strain.L2norm();

  if (strain_norm != 0.0)
    return strain_v / strain_norm;
  else
    return -std::sqrt(3.0);
}

RankTwoTensor
LynxDamageDeformation::rotatedElasticStrain(const RankTwoTensor & elastic_strain)
{
  const Real strain_v = elastic_strain.trace();
  RankTwoTensor strain_rot = spinRotation(elastic_strain.deviatoric());
  strain_rot.addIa(strain_v / 3.0);

  return strain_rot;
}
