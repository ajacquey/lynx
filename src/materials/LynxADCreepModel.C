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

#include "LynxADCreepModel.h"

registerADMooseObject("LynxApp", LynxADCreepModel);

defineADValidParams(
    LynxADCreepModel,
    LynxADMaterialBase,
    params.addClassDescription("Base class for the viscous correction of a visco-elastic rheology.");
    params.set<bool>("compute") = false;
    params.suppressParameter<bool>("compute");
    params.addCoupledVar("temperature", "The temperature variable.");
    // Creep parameters
    params.addParam<std::vector<Real>>("A_diffusion",
                                       "List of pre-exponential factors for diffusion creep.");
    params.addParam<std::vector<Real>>("E_diffusion",
                                       "List of activation energy for diffusion creep.");
    params.addParam<std::vector<Real>>("V_diffusion",
                                       "List of activation (molar) volume for diffusion creep.");
    params.addParam<std::vector<Real>>("A_dislocation",
                                       "List of pre-exponential factors for dislocation creep.");
    params.addParam<std::vector<Real>>("E_dislocation",
                                       "List of activation energy for dislocation creep.");
    params.addParam<std::vector<Real>>("V_dislocation",
                                       "List of activation (molar) volume for dislocation creep.");
    params.addParam<std::vector<Real>>("n_dislocation",
                                       "List of power law exponents for dislocation creep.");
    params.addParam<Real>("gas_constant", 8.3144621, "The universal gas constant");
    params.addParam<std::vector<Real>>(
        "eta_min",
        "The lower threshold for the effective viscosity for purely viscous deformation.");
    params.addParam<std::vector<Real>>(
        "eta_max",
        "The upper threshold for the effective viscosity for purely viscous deformation.");
    MooseEnum viscous_upt("newton=0 brent=1 newton_safe=2", "brent");
    params.addParam<MooseEnum>(
        "viscous_update",
        viscous_upt,
        "The type of update for the viscous strain (Newton, Brent or Newton-safe method).");
    params.addParam<Real>("tolerance", 1.0e-10, "Tolerance for the root finding algorithms.");
    params.addRangeCheckedParam<unsigned int>(
        "max_iterations",
        100,
        "max_iterations>0",
        "Maximum number of iterations in the root finding algorithms."););

template <ComputeStage compute_stage>
LynxADCreepModel<compute_stage>::LynxADCreepModel(const InputParameters & parameters)
  : LynxADMaterialBase<compute_stage>(parameters),
    _coupled_temp(isCoupled("temperature")),
    _temp(_coupled_temp ? adCoupledValue("temperature") : adZeroValue()),
    // Creep parameters
    _has_diffusion_creep(isParamValid("A_diffusion")),
    _A_diffusion(_has_diffusion_creep ? this->getLynxParam("A_diffusion")
                                      : std::vector<Real>(_n_composition, 0.0)),
    _E_diffusion((_has_diffusion_creep && isParamValid("E_diffusion"))
                     ? this->getLynxParam("E_diffusion")
                     : std::vector<Real>(_n_composition, 0.0)),
    _V_diffusion((_has_diffusion_creep && isParamValid("V_diffusion"))
                     ? this->getLynxParam("V_diffusion")
                     : std::vector<Real>(_n_composition, 0.0)),
    _has_dislocation_creep(isParamValid("A_dislocation")),
    _A_dislocation(_has_dislocation_creep ? this->getLynxParam("A_dislocation")
                                          : std::vector<Real>(_n_composition, 0.0)),
    _n_dislocation((_has_dislocation_creep && isParamValid("n_dislocation"))
                       ? this->getLynxParam("n_dislocation")
                       : std::vector<Real>(_n_composition, 1.0)),
    _E_dislocation((_has_dislocation_creep && isParamValid("E_dislocation"))
                       ? this->getLynxParam("E_dislocation")
                       : std::vector<Real>(_n_composition, 0.0)),
    _V_dislocation((_has_dislocation_creep && isParamValid("V_dislocation"))
                       ? this->getLynxParam("V_dislocation")
                       : std::vector<Real>(_n_composition, 0.0)),
    _gas_constant(getParam<Real>("gas_constant")),
    _eta_min(isParamValid("eta_min") ? this->getLynxParam("eta_min")
                                     : std::vector<Real>(_n_composition, 0.0)),
    _eta_max(isParamValid("eta_max") ? this->getLynxParam("eta_max")
                                     : std::vector<Real>(_n_composition, 1.0e+99)),
    _viscous_update(getParam<MooseEnum>("viscous_update")),
    _tol(getParam<Real>("tolerance")),
    _itmax(getParam<unsigned int>("max_iterations")),
    // Creep properties
    _eta_eff(declareADProperty<Real>("effective_viscosity")),
    _viscous_strain_incr(declareADProperty<RankTwoTensor>("viscous_strain_increment"))
{
}

template <ComputeStage compute_stage>
void
LynxADCreepModel<compute_stage>::setQp(unsigned int qp)
{
  _qp = qp;
}

template <ComputeStage compute_stage>
void
LynxADCreepModel<compute_stage>::creepUpdate(ADRankTwoTensor & stress_dev,
                                                 const ADReal & pressure,
                                                 const ADReal & G,
                                                 ADRankTwoTensor & elastic_strain_incr)
{
  const ADReal tau_II = std::sqrt(0.5) * stress_dev.L2norm();
  const ADRankTwoTensor flow_dir = (tau_II != 0.0) ? stress_dev / tau_II : ADRankTwoTensor();

  ADReal delta_e_II = viscousIncrement(pressure, tau_II, G);

  _viscous_strain_incr[_qp] = delta_e_II * flow_dir;
  stress_dev -= 2.0 * G * delta_e_II * flow_dir;
  elastic_strain_incr -= _viscous_strain_incr[_qp];
}

template <ComputeStage compute_stage>
ADReal
LynxADCreepModel<compute_stage>::viscousIncrement(const ADReal & pressure, const ADReal & tau_II, const ADReal & G)
{
  if (tau_II == 0.0)
  {
    _eta_eff[_qp] = 0.0;
    return 0.0;
  }

  // Map creep parameters
  initCreepParameters(pressure);

  if (_has_diffusion_creep && !_has_dislocation_creep) // DIRECT update
  {
    _eta_eff[_qp] = 1.0 / (2.0 * _A_diff);
    _eta_eff[_qp] = std::min(std::max(_eta_eff[_qp], this->averageProperty(_eta_min)), this->averageProperty(_eta_max));
    
    ADReal maxwell_time = _eta_eff[_qp] / G;

    return _dt / maxwell_time / (1.0 + _dt / maxwell_time) * tau_II / (2.0 * G);
  }
  else // ITERATIVE update
  {
    _tau_II_tr = tau_II;
    _G = G;

    ADReal delta_e_II = 0.0;
    switch (_viscous_update)
    {
      case 0: // NEWTON
        delta_e_II = newtonRoot();
        break;
      case 1: // BRENT
      {
        ADReal xneg = _tau_II_tr / (2.0 * _G);
        ADReal xpos = 0.0;
        delta_e_II = brentRoot(xneg, xpos);
        break;
      }
      case 2: // NEWTON SAFE
      {
        ADReal xneg = _tau_II_tr / (2.0 * _G);
        ADReal xpos = 0.0;
        delta_e_II = safeNewtonRoot(xneg, xpos);
        break;
      }
      default:
        mooseError("LynxCreepModel: unknown viscous_update method. Available ones are 'newton', 'brent' or 'newton_safe'!");
    }

    _eta_eff[_qp] = (delta_e_II != 0.0) ? (_tau_II_tr - 2.0 * _G * delta_e_II) * _dt / (2.0 * delta_e_II) : 0.0;

    if ((_eta_eff[_qp] < this->averageProperty(_eta_min)) || (_eta_eff[_qp] > this->averageProperty(_eta_max)))
    {
      _eta_eff[_qp] = std::min(std::max(_eta_eff[_qp], this->averageProperty(_eta_min)), this->averageProperty(_eta_max));
      ADReal maxwell_time = _eta_eff[_qp] / G;
      delta_e_II = _dt / maxwell_time / (1.0 + _dt / maxwell_time) * tau_II / (2.0 * G);
    }
    
    return delta_e_II;
  }  
}

template <ComputeStage compute_stage>
void
LynxADCreepModel<compute_stage>::initCreepParameters(const ADReal & pressure)
{
  _A_diff = this->averageProperty(_A_diffusion);
  _E_diff = this->averageProperty(_E_diffusion);
  _V_diff = this->averageProperty(_V_diffusion);

  _A_disl = this->averageProperty(_A_dislocation);
  _n_disl = this->averageProperty(_n_dislocation);
  _E_disl = this->averageProperty(_E_dislocation);
  _V_disl = this->averageProperty(_V_dislocation);

  if (_coupled_temp)
  { 
    ADReal RT = _gas_constant * _temp[_qp];
    _A_diff *= std::exp(-(_E_diff + pressure * _V_diff) / RT);
    _A_disl *= std::exp(-(_E_disl + pressure * _V_disl) / RT);
  }
}

template <ComputeStage compute_stage>
ADReal
LynxADCreepModel<compute_stage>::iterativeResidual(const ADReal & x)
{
  return _A_diff * (_tau_II_tr - 2.0 * _G * x) * _dt + _A_disl * std::pow(_tau_II_tr - 2.0 * _G * x, _n_disl) * _dt - x;
}

template <ComputeStage compute_stage>
ADReal
LynxADCreepModel<compute_stage>::iterativeResidualDerivative(const ADReal & x)
{
  return -2.0 * _G * _A_diff * _dt - 2.0 * _n_disl * _G * _A_disl * std::pow(_tau_II_tr - 2.0 * _G * x, _n_disl - 1.0) * _dt - 1.0;
}


template <ComputeStage compute_stage>
ADReal
LynxADCreepModel<compute_stage>::newtonRoot()
{
  mooseError("LynxCreepModel: 'newton' update for viscous strain is not yet implemented. Please use 'brent' or 'safe_newton' instead.");
  return 0.0;
}

// Finding root of function using the Van Wijngaarden-Dekker-Brent method or Brent method in short.
// The values x1 and x2 must bracket the root (f(x1) * f(x2) < 0)
// See 9.3 of Numerical Recipes. The Art of Scientific Computing, 3rd Edition, 2007, (C++ code)
template <ComputeStage compute_stage>
ADReal
LynxADCreepModel<compute_stage>::brentRoot(const ADReal & x1, const ADReal x2)
{
  // Machine floating-point precision
  const Real eps = std::numeric_limits<Real>::epsilon();
  // Declare some stuff here
  ADReal a = x1, b = x2, c = x2, d, e, fa = iterativeResidual(a), fb = iterativeResidual(b), fc = fb, p,
       q, r, s, tol1, xm;

  // Chek if x1 and x2 bracket a root
  if (fa * fb > 0.0) // fa and fb are of the same sign
    throw MooseException("LynxCreepModel: in Brent's method, the points x1 and x2 must "
                         "bracket a root of the function!\n");
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
      ADReal min1 = 3.0 * xm * q - std::abs(tol1 * q);
      ADReal min2 = std::abs(e * q);
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
    fb = iterativeResidual(b);
  }
  throw MooseException(
      "LynxCreepModel: maximum number of iterations exceeded in 'brentRoot'!");
}

// Finding root of function using a safe Newton-Raphson method (combinsaison of Newton-Raphson and
// bisection). The values x1 and x2 must bracket the root (f(x1) * f(x2) < 0) See 9.4 of Numerical
// Recipes. The Art of Scientific Computing, 3rd Edition, 2007, (C++ code)
template <ComputeStage compute_stage>
ADReal
LynxADCreepModel<compute_stage>::safeNewtonRoot(const ADReal & x1, const ADReal x2)
{
  // Machine floating-point precision
  const Real eps = std::numeric_limits<Real>::epsilon();
  // Declare some stuff here
  ADReal xh, xl;
  ADReal fl = iterativeResidual(x1);
  ADReal fh = iterativeResidual(x2);
  if (fl * fh > 0.0) // fl and fh are of the same sign
    throw MooseException("LynxCreepModel: in 'safeNewtonRoot', the points x1 and x2 "
                         "must bracket a root of the function!\n");

  if (fl == 0.0)
    return x1;
  if (fh == 0.0)
    return x2;
  // Orient the search so that f(x1) < 0
  if (fl < 0.0)
  {
    xl = x1;
    xh = x2;
  }
  else
  {
    xh = x1;
    xl = x2;
  }
  // Initialize the guess for root
  ADReal rts = 0.5 * (x1 + x2);
  // The "stepsize before last"
  ADReal dxold = std::abs(x2 - x1);
  // and the last step
  ADReal dx = dxold;
  ADReal f = iterativeResidual(rts);
  ADReal df = iterativeResidualDerivative(rts);
  // Loop over allowed iterations
  for (unsigned int iter = 0; iter < _itmax; ++iter)
  {
    // Bisect if Newton is out of range or not decreasing fast enough
    if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) ||
        (std::abs(2.0 * f) > std::abs(dxold * df)))
    {
      dxold = dx;
      dx = 0.5 * (xh - xl);
      rts = xl + dx;
      // Change in root is negligible.
      if (xl == rts)
        return rts;
    }
    else
    {
      // Newton step acceptable. Take it.
      dxold = dx;
      dx = f / df;
      ADReal temp = rts;
      rts -= dx;
      // Change in root is negligible.
      if (temp == rts)
        return rts;
    }
    // Convergence criterion
    if (std::abs(dx) < (2.0 * eps * std::abs(rts) + 0.5 * _tol))
      return rts;
    f = iterativeResidual(rts);
    df = iterativeResidualDerivative(rts);
    // The one new function evaluation per iteration.
    // Maintain the bracket on the root
    if (f < 0.0)
      xl = rts;
    else
      xh = rts;
  }
  throw MooseException(
      "LynxCreepModel: maximum number of iterations exceeded in 'safeNewtonRoot'!");
}

// adBaseClass(LynxADCreepModel);