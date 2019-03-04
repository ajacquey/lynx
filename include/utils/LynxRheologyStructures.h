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

#ifndef LYNXRHEOLOGYSTRUCTURES_H
#define LYNXRHEOLOGYSTRUCTURES_H

struct diffusion_creep
{
  Real _A0;
  Real _A;
  Real _act_energy;
  Real _act_vol;
  // Real _reference_gs;
  // Real _gs;
  // Real _n_gs;
  // Real _COH;
  // Real _n_COH;
  diffusion_creep() : _A0(0.0), _A(0.0), _act_energy(0.0), _act_vol(0.0)
  // _reference_gs(1.0),
  // _gs(1.0),
  // _n_gs(0.0),
  // _COH(1.0),
  // _n_COH(0.0)
  {
  }
  void fill(const Real A, const Real E, const Real V)
  // const Real gs_ref,
  // const Real gs,
  // const Real n_gs,
  // const Real COH,
  // const Real n_COH)
  {
    _A0 = A;
    _act_energy = E;
    _act_vol = V;
    // _reference_gs = gs_ref;
    // _gs = gs;
    // _n_gs = n_gs;
    // _COH = COH;
    // _n_COH = n_COH;
  }
  Real creepRate(const Real eqv_stress)
  {
    Real s_II = eqv_stress / std::sqrt(3.0);
    return 2.0 / std::sqrt(3.0) * _A * s_II;
  }
  Real creepRateDerivative(const Real /*eqv_stress*/) { return 2.0 / 3.0 * _A; }
  Real oneOnViscosity() { return 2.0 * _A; }
};

struct dislocation_creep
{
  Real _A0;
  Real _A;
  Real _n_exp;
  Real _act_energy;
  Real _act_vol;
  dislocation_creep() : _A0(0.0), _A(0.0), _n_exp(1.0), _act_energy(0.0), _act_vol(0.0) {}
  void fill(const Real A, const Real n, const Real E, const Real V)
  {
    _A0 = A;
    _act_energy = E;
    _act_vol = V;
    _n_exp = n;
  }
  Real creepRate(const Real eqv_stress)
  {
    Real s_II = eqv_stress / std::sqrt(3.0);
    return 2.0 / std::sqrt(3.0) * _A * std::pow(s_II, _n_exp);
  }
  Real creepRateDerivative(const Real eqv_stress)
  {
    Real s_II = eqv_stress / std::sqrt(3.0);
    return 2.0 / 3.0 * _A * _n_exp * std::pow(s_II, _n_exp - 1.0);
  }
  Real oneOnViscosity(const Real strain_rate_II)
  {
    return 2.0 * std::pow(_A, 1.0 / _n_exp) * std::pow(strain_rate_II, 1.0 - 1.0 / _n_exp);
  }
  Real oneOnViscosityDerivative(const Real strain_rate_II)
  {
    return 2.0 * std::pow(_A, 1.0 / _n_exp) * (1.0 - 1.0 / _n_exp) *
           std::pow(strain_rate_II, -1.0 / _n_exp);
  }
};

struct iterative_viscous
{
  diffusion_creep * _diff_creep;
  dislocation_creep * _disl_creep;
  Real _eqv_stress;
  Real _G;
  Real _dt;
  bool _has_L;
  bool _has_N;
  iterative_viscous() : _eqv_stress(0.0), _G(0.0), _dt(1.0), _has_L(false), _has_N(false) {}
  void fill(const Real eqv_stress, const Real G, const Real dt)
  {
    _eqv_stress = eqv_stress;
    _G = G;
    _dt = dt;
    _has_L = false;
    _has_N = false;
  }
  void fillDiff(diffusion_creep * diff_creep)
  {
    _diff_creep = diff_creep;
    _has_L = true;
  }
  void fillDisl(dislocation_creep * disl_creep)
  {
    _disl_creep = disl_creep;
    _has_N = true;
  }
  Real func(const Real x)
  {
    Real eqv_stress = _eqv_stress - 3.0 * _G * x * _dt;
    if (eqv_stress < 1.0e-20)
      return -x;
    Real L_rate = 0.0, N_rate = 0.0;
    if (_has_L)
      L_rate = _diff_creep->creepRate(eqv_stress);
    if (_has_N)
      N_rate = _disl_creep->creepRate(eqv_stress);

    Real creep_rate = L_rate + N_rate;
    return creep_rate - x;
  }
  Real deriv(const Real x)
  {
    Real eqv_stress = _eqv_stress - 3.0 * _G * x * _dt;
    Real L_rate_deriv = 0.0, N_rate_deriv = 0.0;
    if (_has_L)
      L_rate_deriv = _diff_creep->creepRateDerivative(eqv_stress) * 3.0 * _G * _dt;
    if (_has_N)
      N_rate_deriv = _disl_creep->creepRateDerivative(eqv_stress) * 3.0 * _G * _dt;

    Real creep_rate_deriv = L_rate_deriv + N_rate_deriv;
    return creep_rate_deriv - 1.0;
  }
};

struct plasticity
{
  Real _alpha;
  Real _alpha_0;
  Real _alpha_res;
  Real _k;
  Real _k_0;
  Real _k_res;
  Real _beta;
  Real _H;
  Real _eta; // actually one on eta
  Real _intnl;
  Real _intnl_0;
  Real _intnl_lim;
  plasticity()
    : _alpha(0.0),
      _alpha_0(0.0),
      _alpha_res(0.0),
      _k(0.0),
      _k_0(0.0),
      _k_res(0.0),
      _beta(0.0),
      _H(0.0),
      _eta(0.0),
      _intnl(0.0),
      _intnl_0(0.0),
      _intnl_lim(0.0)
  {
  }
  void fill(const Real alpha_0,
            const Real alpha_res,
            const Real k_0,
            const Real k_res,
            const Real beta,
            const Real intnl_0,
            const Real intnl_lim,
            const Real eta)
  {
    _alpha_0 = alpha_0;
    _alpha_res = alpha_res;
    _k_0 = k_0;
    _k_res = k_res;
    _beta = beta;
    _intnl_0 = intnl_0;
    _intnl_lim = intnl_lim;
    _eta = eta;
  }
};

struct damage_plasticity
{
  Real _xi0;
  Real _gamma;
  Real _p_cr;
  Real _eta_p; // actually one on eta
  Real _eta_d; // actually one on eta
  Real _alpha0;
  Real _p_tr;        // trial pressure
  Real _q_tr;        // trial eqv_stress
  Real _p_r;         // reference pressure for convex yield
  Real _q_r;         // reference eqv_stress for convex yield
  Real _rho_tr;      // trial distance to reference point
  Real _rp;          // direction for projection
  Real _rq;          // direction for projection
  Real _xi_cr;       // critical strain ratio for capped yield
  Real _alpha2;      // square of alpha parameter for capped yield
  Real _dxi_cr_dp;   // derivative wrt pressure of the critical strain ratio
  Real _dxi_cr_dq;   // derivative wrt eqv_stress of the critical strain ratio
  Real _dmu2_dxi_cr; // derivative wrt the critical strain ratio of the square of alpha parameter
  damage_plasticity()
    : _xi0(-std::sqrt(3.0)),
      _gamma(0.0),
      _p_cr(0.0),
      _eta_p(0.0),
      _eta_d(0.0),
      _alpha0(0.0),
      _p_tr(0.0),
      _q_tr(0.0),
      _p_r(0.0),
      _q_r(0.0),
      _rho_tr(0.0),
      _rp(0.0),
      _rq(0.0),
      _xi_cr(-std::sqrt(3.0)),
      _alpha2(0.0),
      _dxi_cr_dp(0.0),
      _dxi_cr_dq(0.0),
      _dmu2_dxi_cr(0.0)
  {
  }
  void fill(const Real xi0, const Real gamma, const Real p_cr, const Real eta_p, const Real eta_d)
  {
    _xi0 = xi0;
    _gamma = gamma;
    _p_cr = p_cr;
    _eta_p = eta_p;
    _eta_d = eta_d;
  }
};

#endif // LYNXRHEOLOGYSTRUCTURES_H
