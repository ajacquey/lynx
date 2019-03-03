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
  InputParameters params = validParams<LynxDeformation>();
  params.addClassDescription("Class calculating the deformation of a material following a "
                             "damage viscoelasto-(visco)plastic rheology.");
  // Coupled variables
  params.addCoupledVar("damage", "The damage variable.");
  params.addCoupledVar("porosity", "The porosity variable.");
  // Strain parameters
  // Elastic moduli parameters
  params.addParam<std::vector<Real>>("damage_modulus",
                                     "The third elastic damage modulus for the damage heology.");
  // Creep parameters
  // Plastic parameters
  params.addParam<std::vector<Real>>("critical_pressure",
                                     "The critical pressure for capped plastic models.");
  return params;
}

LynxDamageDeformation::LynxDamageDeformation(const InputParameters & parameters)
  : LynxDeformation(parameters),
    // Coupled variables
    _coupled_dam(isCoupled("damage")),
    _damage(_coupled_dam ? coupledValue("damage") : _zero),
    _damage_old(_coupled_dam ? coupledValueOld("damage") : _zero),
    _coupled_phi(isCoupled("porosity")),
    _porosity(_coupled_phi ? coupledValue("porosity") : _zero),
    // Strain parameters
    // Elastic moduli parameters
    _damage_modulus(isParamValid("damage_modulus") ? getLynxParam<Real>("damage_modulus")
                                                   : std::vector<Real>(_n_composition, 0.0)),
    // Creep parameters
    // Plastic parameters
    _critical_pressure((_has_plasticity && isParamValid("critical_pressure"))
                           ? getLynxParam<Real>("critical_pressure")
                           : std::vector<Real>(_n_composition, 0.0)),
    // Rheology boolean
    // Strain properties
    _elastic_strain(declareProperty<RankTwoTensor>("elastic_strain")),
    _elastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("elastic_strain"))
// Viscous properties
// Plastic properties
// Stress properties
{
}

void
LynxDamageDeformation::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
  LynxDeformation::initQpStatefulProperties();
}

void
LynxDamageDeformation::computeQpDeformation()
{
  // Update elastic moduli
  elasticModuli();

  // Update the volumetric part of the deformation
  Real pressure = -_stress_old[_qp].trace() / 3.0;
  volumetricDeformation(pressure);

  // Update the deviatoric part of the deformation
  RankTwoTensor stress_dev = spinRotation(_stress_old[_qp].deviatoric());
  deviatoricDeformation(pressure, stress_dev);

  // Plastic correction
  if (_has_plasticity && _G[_qp] != 0.0)
    plasticCorrection(pressure, stress_dev);

  // Form the total stress tensor
  _stress[_qp] = stress_dev;
  _stress[_qp].addIa(-pressure);

  // Update tangent operator modulus
  if (_fe_problem.currentlyComputingJacobian())
    tangentOperator();
}

void
LynxDamageDeformation::elasticModuli()
{
  // Elastic strain ratio
  const Real xi_old = strainRatio(_elastic_strain_old[_qp]);

  // Critical strain strain ratio
  _xi0 =

  // Bulk modulus
  _K[_qp] = _has_bulk_modulus ? averageProperty(_bulk_modulus) : 0.0;

  // Shear modulus
  _G[_qp] = _has_shear_modulus ? averageProperty(_shear_modulus) : 0.0;

  // Initialize inelastic increment
  _viscous_strain_incr[_qp].zero();
  _plastic_strain_incr[_qp].zero();

  // Initialize stuff
  _stress_corr_v = 1.0;
  _dq_dq_tr_v = 1.0;
  _stress_corr_p = 1.0;
  _dp_dp_tr_p = 1.0;
  _dp_dq_tr_p = 0.0;
  _dq_dp_tr_p = 0.0;
  _dq_dq_tr_p = 1.0;
  _plastic_yield_function[_qp] = -1.0;
}

void
LynxDamageDeformation::plasticCorrection(Real & pressure, RankTwoTensor & stress_dev)
{
  // Check if yield is capped
  _p_cr = averageProperty(_critical_pressure);
  if (_p_cr > 0.0)
    convexPlasticCorrection(pressure, stress_dev);
  else if (_p_cr < 0.0)
    mooseError("LynxDamageDeformation: 'critical_pressure' must be positive!");
  else
  {
    Real eqv_stress = std::sqrt(1.5) * stress_dev.L2norm();
    RankTwoTensor flow_dir = (eqv_stress != 0.0) ? stress_dev / eqv_stress : RankTwoTensor();

    // Plastic correction
    Real plastic_incr = plasticIncrement(pressure, eqv_stress);

    _plastic_strain_incr[_qp] = 1.5 * plastic_incr * flow_dir;
    _plastic_strain_incr[_qp].addIa(_plasticity->_beta * plastic_incr);
    stress_dev -= 3.0 * _G[_qp] * plastic_incr * flow_dir;
    pressure += _K[_qp] * _plasticity->_beta * plastic_incr;
  }
}

void
LynxDamageDeformation::convexPlasticCorrection(Real & pressure, RankTwoTensor & stress_dev)
{
  RankTwoTensor stress_dev_tr = stress_dev;
  Real eqv_stress = std::sqrt(1.5) * stress_dev.L2norm();
  RankTwoTensor flow_dir = (eqv_stress != 0.0) ? stress_dev / eqv_stress : RankTwoTensor();

  Real sin_phi = std::sin(averageProperty(_friction_angle_0) * libMesh::pi / 180.0);
  _xi0 = -std::sqrt(3.0) / std::sqrt(1.0 + 1.5 * Utility::pow<2>(_K[_qp] / _G[_qp] * sin_phi));
  Real k = std::sqrt(3.0) * averageProperty(_cohesion_0) *
           std::cos(averageProperty(_friction_angle_0) * libMesh::pi / 180.0);
  _p_k = -k / std::sqrt(2.0 * Utility::pow<2>(_G[_qp] / _K[_qp]) * (3.0 - Utility::pow<2>(_xi0)) /
                        Utility::pow<2>(_xi0));
  Real yield_squared = convexPlasticYield2(pressure, eqv_stress);

  if (yield_squared > 0.0)
  {
    Real vp_correction = 1.0;
    Real one_on_eta = averageProperty(_one_on_plastic_eta);

    // Trial state
    _p_tr = pressure;
    _q_tr = eqv_stress;

    // Reference point
    _p_r = convexReferencePressure(_p_tr, _q_tr);
    _q_r = 0.0;

    _rho_tr = std::sqrt(Utility::pow<2>(_p_tr - _p_r) + Utility::pow<2>(_q_tr - _q_r));
    _rp = (_p_tr - _p_r) / _rho_tr;
    _rq = (_q_tr - _q_r) / _rho_tr;

    // Find projection on convex yield
    Real rho_neg = 0.0;
    Real rho_pos = _rho_tr;
    Real rho_0 = getConvexProjection(rho_neg, rho_pos);

    // Update projections
    Real p_0 = _p_r + rho_0 * _rp;
    Real q_0 = _q_r + rho_0 * _rq;

    // Plastic strain
    Real vol_strain_p = (_p_tr - p_0) / _K[_qp];
    Real eqv_strain_p = (_q_tr - q_0) / (3.0 * _G[_qp]);

    // Viscoplastic correction
    if (one_on_eta != 0.0)
    {
      Real F_tr = _rho_tr - rho_0;
      Real delta_p = std::sqrt(Utility::pow<2>(vol_strain_p) + Utility::pow<2>(eqv_strain_p));

      vp_correction = (F_tr * one_on_eta * _dt / delta_p) / (1.0 + F_tr * one_on_eta * _dt / delta_p);

      vol_strain_p *= vp_correction;
      eqv_strain_p *= vp_correction;
    }

    _plastic_strain_incr[_qp] = 1.5 * eqv_strain_p * flow_dir;
    _plastic_strain_incr[_qp].addIa(-vol_strain_p / 3.0);
    pressure -= _K[_qp] * vol_strain_p;
    stress_dev -= 3.0 * _G[_qp] * eqv_strain_p * flow_dir;
    eqv_stress = std::sqrt(1.5) * stress_dev.L2norm();
    _plastic_yield_function[_qp] =
        std::sqrt(Utility::pow<2>(pressure - _p_r) + Utility::pow<2>(eqv_stress - _q_r)) - rho_0;

    // Tangent operator
    _stress_corr_p = (eqv_stress != 0.0) ? eqv_stress / _q_tr : 1.0;
    Real dF2_dp0 = dConvexPlasticYield2_dp(p_0, q_0) / dConvexPlasticYield2(rho_0);
    Real dF2_dq0 = dConvexPlasticYield2_dq(p_0, q_0) / dConvexPlasticYield2(rho_0);

    _dp_dp_tr_p =
        _dp_r_dp_tr + (1.0 - vp_correction) + vp_correction * rho_0 / _rho_tr * dF2_dq0 * _rq;
    _dp_dq_tr_p = _dp_r_dq_tr - vp_correction * rho_0 / _rho_tr * dF2_dq0 * _rp;
    _dq_dp_tr_p = -vp_correction * rho_0  / _rho_tr * dF2_dp0 * _rq;
    _dq_dq_tr_p = (1.0 - vp_correction) + vp_correction * rho_0  / _rho_tr * dF2_dp0 * _rp;
  }
}

Real
LynxDamageDeformation::convexReferencePressure(const Real & p_tr, const Real & q_tr)
{
  Real mu0 = std::sqrt(2.0 * Utility::pow<2>(_G[_qp] / _K[_qp]) * (3.0 - Utility::pow<2>(_xi0)) /
                       Utility::pow<2>(_xi0));
  Real k = std::sqrt(3.0) * averageProperty(_cohesion_0) *
           std::cos(averageProperty(_friction_angle_0) * libMesh::pi / 180.0);

  Real p_proj =
      (1.0 - Utility::pow<2>(mu0) * _K[_qp] / (3.0 * _G[_qp] + Utility::pow<2>(mu0) * _K[_qp])) *
          p_tr +
      mu0 * _K[_qp] * (q_tr - k) / (3.0 * _G[_qp] + Utility::pow<2>(mu0) * _K[_qp]);

  // _dp_r_dp_tr = ((p_proj > 0.0) && (p_proj < _p_cr)) ? (1.0 - Utility::pow<2>(mu0) * _K[_qp] /
  //                                            (3.0 * _G[_qp] + Utility::pow<2>(mu0) * _K[_qp]))
  //                               : 0.0;
  // _dp_r_dq_tr = ((p_proj > 0.0) && (p_proj < _p_cr)) ? mu0 * _K[_qp] / (3.0 * _G[_qp] +
  // Utility::pow<2>(mu0) * _K[_qp]) : 0.0;
  return ((p_proj > 0.0) && (p_proj < _p_cr)) ? p_proj : 2.0 * _p_cr / 3.0;
}

Real
LynxDamageDeformation::convexPlasticYield2(const Real & rho)
{
  Real p = _p_r + rho * _rp;
  Real q = _q_r + rho * _rq;

  return convexPlasticYield2(p, q);
}

Real
LynxDamageDeformation::dConvexPlasticYield2(const Real & rho)
{
  Real p = _p_r + rho * _rp;
  Real q = _q_r + rho * _rq;

  return dConvexPlasticYield2_dp(p, q) * _rp + dConvexPlasticYield2_dq(p, q) * _rq;
}

Real
LynxDamageDeformation::convexPlasticYield2(const Real & pressure, const Real & eqv_stress)
{
  if ((pressure == 0.0) && (eqv_stress == 0.0))
    return 0.0;

  Real pq_term =
      (pressure != 0.0)
          ? 1.0 / (1.0 + 0.5 * Utility::pow<2>(_K[_qp] * eqv_stress / (_G[_qp] * pressure)))
          : 0.0;
  Real xi_cr =
      (pressure > 0.0) ? _xi0 - (std::sqrt(3.0) + _xi0) * pressure / _p_cr * pq_term : _xi0;

  Real mu2 = (xi_cr != 0.0) ? 2.0 * Utility::pow<2>(_G[_qp] / _K[_qp]) *
                                  (3.0 - Utility::pow<2>(xi_cr)) / Utility::pow<2>(xi_cr)
                            : 0.0;

  return Utility::pow<2>(eqv_stress) - mu2 * Utility::pow<2>(pressure - _p_k) *
                                           ((0.0 < pressure - _p_k) - (pressure - _p_k < 0.0));
}

Real
LynxDamageDeformation::dConvexPlasticYield2_dp(const Real & pressure, const Real & eqv_stress)
{
  if ((pressure == 0.0) && (eqv_stress == 0.0))
    return 0.0;

  Real pq_term =
      (pressure != 0.0)
          ? 1.0 / (1.0 + 0.5 * Utility::pow<2>(_K[_qp] * eqv_stress / (_G[_qp] * pressure)))
          : 0.0;
  Real xi_cr =
      (pressure > 0.0) ? _xi0 - (std::sqrt(3.0) + _xi0) * pressure / _p_cr * pq_term : _xi0;
  Real dxi_cr_dp = (pressure != 0.0) ? -(_xi0 - xi_cr) / pressure * (3.0 - 2.0 * pq_term) : 0.0;

  Real mu2 = (xi_cr != 0.0) ? 2.0 * Utility::pow<2>(_G[_qp] / _K[_qp]) *
                                  (3.0 - Utility::pow<2>(xi_cr)) / Utility::pow<2>(xi_cr)
                            : 0.0;
  Real dmu2_dxi_cr =
      (xi_cr != 0.0) ? -12.0 * Utility::pow<2>(_G[_qp] / _K[_qp]) / Utility::pow<3>(xi_cr) : 0.0;

  return -(dmu2_dxi_cr * dxi_cr_dp * (pressure - _p_k) + 2.0 * mu2) * (pressure - _p_k) *
         ((0.0 < pressure - _p_k) - (pressure - _p_k < 0.0));
}

Real
LynxDamageDeformation::dConvexPlasticYield2_dq(const Real & pressure, const Real & eqv_stress)
{
  if ((pressure == 0.0) && (eqv_stress == 0.0))
    return 0.0;

  Real pq_term =
      (pressure != 0.0)
          ? 1.0 / (1.0 + 0.5 * Utility::pow<2>(_K[_qp] * eqv_stress / (_G[_qp] * pressure)))
          : 0.0;
  Real xi_cr =
      (pressure > 0.0) ? _xi0 - (std::sqrt(3.0) + _xi0) * pressure / _p_cr * pq_term : _xi0;
  Real dxi_cr_dq = (eqv_stress != 0.0) ? 2.0 * (_xi0 - xi_cr) / eqv_stress * (1.0 - pq_term) : 0.0;
  Real dmu2_dxi_cr =
      (xi_cr != 0.0) ? -12.0 * Utility::pow<2>(_G[_qp] / _K[_qp]) / Utility::pow<3>(xi_cr) : 0.0;

  return 2.0 * eqv_stress -
         dmu2_dxi_cr * dxi_cr_dq *
             Utility::pow<2>(pressure -
                             _p_k); // *
                                    // ((0.0 < pressure - _p_k) - (pressure - _p_k < 0.0));
}

Real
LynxDamageDeformation::getConvexProjection(const Real & x1, const Real & x2)
{
  // Machine floating-point precision
  const Real eps = std::numeric_limits<Real>::epsilon();
  // Declare some stuff here
  Real a = x1, b = x2, c = x2, d, e, fa = convexPlasticYield2(a), fb = convexPlasticYield2(b), fc,
       p, q, r, s, tol1, xm;
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
    fb = convexPlasticYield2(b);
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
