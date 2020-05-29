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

registerMooseObject("LynxApp", LynxADCreepModel);

InputParameters
LynxADCreepModel::validParams()
{
  InputParameters params = LynxADMaterialBase::validParams();
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
      "initial_viscosity", "The vector of initial viscosity for purely viscous deformation.");
  params.addParam<Real>("background_strain_rate",
                        "The background strain rate for purely viscous deformation.");
  params.addParam<std::vector<Real>>(
      "eta_min", "The lower threshold for the effective viscosity for purely viscous deformation.");
  params.addParam<std::vector<Real>>(
      "eta_max", "The upper threshold for the effective viscosity for purely viscous deformation.");
  MooseEnum viscous_upt("direct=0 brent=1 newton_safe=2", "direct");
  params.addParam<MooseEnum>(
      "viscous_update",
      viscous_upt,
      "The type of update for the viscous strain (Direct, Brent or Newton-safe method).");
  params.addParam<Real>("tolerance", 1.0e-10, "Tolerance for the root finding algorithms.");
  params.addRangeCheckedParam<unsigned int>(
      "max_iterations",
      100,
      "max_iterations>0",
      "Maximum number of iterations in the root finding algorithms.");
  return params;
}

LynxADCreepModel::LynxADCreepModel(const InputParameters & parameters)
  : LynxADMaterialBase(parameters),
    _coupled_temp(isCoupled("temperature")),
    _temp(_coupled_temp ? adCoupledValue("temperature") : adZeroValue()),
    // Creep parameters
    _has_diffusion_creep(isParamValid("A_diffusion")),
    _A_diffusion(_has_diffusion_creep ? getLynxParam<Real>("A_diffusion")
                                      : std::vector<Real>(_n_composition, 0.0)),
    _E_diffusion((_has_diffusion_creep && isParamValid("E_diffusion"))
                     ? getLynxParam<Real>("E_diffusion")
                     : std::vector<Real>(_n_composition, 0.0)),
    _V_diffusion((_has_diffusion_creep && isParamValid("V_diffusion"))
                     ? getLynxParam<Real>("V_diffusion")
                     : std::vector<Real>(_n_composition, 0.0)),
    _has_dislocation_creep(isParamValid("A_dislocation")),
    _A_dislocation(_has_dislocation_creep ? getLynxParam<Real>("A_dislocation")
                                          : std::vector<Real>(_n_composition, 0.0)),
    _n_dislocation((_has_dislocation_creep && isParamValid("n_dislocation"))
                       ? getLynxParam<Real>("n_dislocation")
                       : std::vector<Real>(_n_composition, 1.0)),
    _E_dislocation((_has_dislocation_creep && isParamValid("E_dislocation"))
                       ? getLynxParam<Real>("E_dislocation")
                       : std::vector<Real>(_n_composition, 0.0)),
    _V_dislocation((_has_dislocation_creep && isParamValid("V_dislocation"))
                       ? getLynxParam<Real>("V_dislocation")
                       : std::vector<Real>(_n_composition, 0.0)),
    _gas_constant(getParam<Real>("gas_constant")),
    _has_background_strain_rate(isParamValid("background_strain_rate")),
    _has_initial_viscosity(_has_background_strain_rate ? false
                                                       : isParamValid("_has_initial_viscosity")),
    _initial_viscosity(_has_initial_viscosity ? getLynxParam<Real>("initial_viscosity")
                                              : std::vector<Real>(_n_composition, 0.0)),
    _background_strain_rate(_has_background_strain_rate ? getParam<Real>("background_strain_rate")
                                                        : 0.0),
    _eta_min(isParamValid("eta_min") ? getLynxParam<Real>("eta_min")
                                     : std::vector<Real>(_n_composition, 0.0)),
    _eta_max(isParamValid("eta_max") ? getLynxParam<Real>("eta_max")
                                     : std::vector<Real>(_n_composition, 1.0e+99)),
    _viscous_update(getParam<MooseEnum>("viscous_update")),
    _tol(getParam<Real>("tolerance")),
    _itmax(getParam<unsigned int>("max_iterations")),
    // Creep properties
    _eta_eff(declareADProperty<Real>("effective_viscosity")),
    _viscous_strain_incr(declareADProperty<RankTwoTensor>("viscous_strain_increment"))
{
}

void
LynxADCreepModel::setQp(unsigned int qp)
{
  _qp = qp;
}

void
LynxADCreepModel::creepUpdate(ADRankTwoTensor & stress_dev,
                              const ADReal & pressure,
                              const Real & G,
                              ADRankTwoTensor & elastic_strain_incr)
{
  if (G != 0.0) // Visco elastic update
  {
    const ADReal tau_II_tr = std::sqrt(0.5) * stress_dev.L2norm();
    const ADRankTwoTensor flow_dir =
        (tau_II_tr != 0.0) ? stress_dev / tau_II_tr : ADRankTwoTensor();

    ADReal delta_e_II = viscousIncrement(pressure, tau_II_tr, G);

    _viscous_strain_incr[_qp] = delta_e_II * flow_dir;
    stress_dev -= 2.0 * G * delta_e_II * flow_dir;
    elastic_strain_incr -= _viscous_strain_incr[_qp];
  }
  else // Stoke
  {
    _viscous_strain_incr[_qp] = elastic_strain_incr;
    // Compute effective viscosity
    _eta_eff[_qp] = computeQpEffectiveViscosity(pressure);

    stress_dev = 2.0 * _eta_eff[_qp] * _viscous_strain_incr[_qp].deviatoric() / _dt;
  }
}

ADReal
LynxADCreepModel::viscousIncrement(const ADReal & pressure,
                                   const ADReal & tau_II_tr,
                                   const Real & G)
{
  if (tau_II_tr == 0.0)
  {
    _eta_eff[_qp] = 0.0;
    return 0.0;
  }

  // Map creep parameters
  initCreepParameters(pressure);

  _tau_II_tr = tau_II_tr;
  _G = G;

  if (_viscous_update == 0) // DIRECT
  {
    _eta_eff[_qp] =
        0.5 * _tau_II_tr / (_A_diff * _tau_II_tr + _A_disl * std::pow(_tau_II_tr, _n_disl));
  }
  else // ITERATIVE
  {
    ADReal tau_II = 0.0;

    switch (_viscous_update)
    {
      case 1: // BRENT
      {
        ADReal xneg = _tau_II_tr;
        ADReal xpos = 0.0;
        tau_II = brentRoot(xneg, xpos);
        break;
      }
      case 2: // NEWTON SAFE
      {
        ADReal xneg = _tau_II_tr;
        ADReal xpos = 0.0;
        tau_II = safeNewtonRoot(xneg, xpos);
        break;
      }
      default:
        mooseError("LynxCreepModel: unknown viscous_update method. Available ones are 'direct', "
                   "'brent' or 'newton_safe'!");
    }
    _eta_eff[_qp] = (tau_II != _tau_II_tr) ? _G * _dt * tau_II / (_tau_II_tr - tau_II) : 0.0;
  }

  if ((_eta_eff[_qp] < averageProperty(_eta_min)) || (_eta_eff[_qp] > averageProperty(_eta_max)))
    _eta_eff[_qp] =
        std::min(std::max(_eta_eff[_qp], averageProperty(_eta_min)), averageProperty(_eta_max));
  ADReal maxwell_time = _eta_eff[_qp] / _G;
  ADReal delta_e_II = (maxwell_time != 0.0) ? _dt / maxwell_time / (1.0 + _dt / maxwell_time) *
                                                  _tau_II_tr / (2.0 * _G)
                                            : 0.0;

  return delta_e_II;
}

ADReal
LynxADCreepModel::computeQpEffectiveViscosity(const ADReal & pressure)
{
  // Map creep parameters
  initCreepParameters(pressure);

  ADReal strain_rate_II = std::sqrt(0.5) * _viscous_strain_incr[_qp].deviatoric().L2norm() / _dt;

  if (_t_step <= 1 && strain_rate_II == 0.0 && _has_initial_viscosity)
    return std::min(std::max(averageProperty(_initial_viscosity), averageProperty(_eta_min)),
                    averageProperty(_eta_max));

  strain_rate_II = (_has_background_strain_rate && _t_step <= 1 && strain_rate_II == 0.0)
                       ? _background_strain_rate
                       : strain_rate_II;

  ADReal one_on_eta_diff = 0.0, one_on_eta_disl = 0.0;

  if (_has_diffusion_creep)
    one_on_eta_diff = computeQpOneOnDiffViscosity(_A_diff);
  if (_has_dislocation_creep)
    one_on_eta_disl = computeQpOneOnDislViscosity(_A_disl, _n_disl, strain_rate_II);

  ADReal eta = 1.0 / (one_on_eta_diff + one_on_eta_disl);

  return std::min(std::max(eta, averageProperty(_eta_min)), averageProperty(_eta_max));
}

ADReal
LynxADCreepModel::computeQpOneOnDiffViscosity(const ADReal A)
{
  return 2.0 * A;
}

ADReal
LynxADCreepModel::computeQpOneOnDislViscosity(const ADReal A, const ADReal n, const ADReal eII)
{
  if ((eII == 0.0) && (n == 1.0))
    return 2.0;
  else
    return 2.0 * std::pow(A, 1.0 / n) * std::pow(eII, 1.0 - 1.0 / n);
}

void
LynxADCreepModel::initCreepParameters(const ADReal & pressure)
{
  _A_diff = averageProperty(_A_diffusion);
  _E_diff = averageProperty(_E_diffusion);
  _V_diff = averageProperty(_V_diffusion);

  _A_disl = averageProperty(_A_dislocation);
  _n_disl = averageProperty(_n_dislocation);
  _E_disl = averageProperty(_E_dislocation);
  _V_disl = averageProperty(_V_dislocation);

  if (_coupled_temp)
  {
    ADReal RT = _gas_constant * _temp[_qp];
    _A_diff *= std::exp(-(_E_diff + pressure * _V_diff) / RT);
    _A_disl *= std::exp(-(_E_disl + pressure * _V_disl) / RT);
  }
}

ADReal
LynxADCreepModel::iterativeResidual(const ADReal & x)
{
  return 2.0 * _G * (_A_diff * x + _A_disl * std::pow(x, _n_disl)) * _dt - (_tau_II_tr - x);
}

ADReal
LynxADCreepModel::iterativeResidualDerivative(const ADReal & x)
{
  return 2.0 * _G * (_A_diff + _n_disl * _A_disl * std::pow(x, _n_disl - 1.0)) * _dt + 1.0;
}

ADReal
LynxADCreepModel::newtonRoot()
{
  mooseError("LynxCreepModel: 'newton' update for viscous strain is not yet implemented. Please "
             "use 'brent' or 'safe_newton' instead.");
  return 0.0;
}

// Finding root of function using the Van Wijngaarden-Dekker-Brent method or Brent method in
// short. The values x1 and x2 must bracket the root (f(x1) * f(x2) < 0) See 9.3 of Numerical
// Recipes. The Art of Scientific Computing, 3rd Edition, 2007, (C++ code)
ADReal
LynxADCreepModel::brentRoot(const ADReal & x1, const ADReal x2)
{
  // Machine floating-point precision
  const Real eps = std::numeric_limits<Real>::epsilon();
  // Declare some stuff here
  ADReal a = x1, b = x2, c = x2, d, e, fa = iterativeResidual(a), fb = iterativeResidual(b),
         fc = fb, p, q, r, s, tol1, xm;

  // Chek if x1 and x2 bracket a root
  if (fa * fb > 0.0) // fa and fb are of the same sign
    throw MooseException("LynxCreepModel: in Brent's method, the points x1 and x2 must "
                         "bracket a root of the function!\nx1 = ",
                         x1,
                         " and f(x1) = ",
                         fa,
                         "\nx2 = ",
                         x2,
                         " and f(x2) = ",
                         fb,
                         "\n");
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
  throw MooseException("LynxCreepModel: maximum number of iterations exceeded in 'brentRoot'!");
}

// Finding root of function using a safe Newton-Raphson method (combinsaison of Newton-Raphson and
// bisection). The values x1 and x2 must bracket the root (f(x1) * f(x2) < 0) See 9.4 of Numerical
// Recipes. The Art of Scientific Computing, 3rd Edition, 2007, (C++ code)
ADReal
LynxADCreepModel::safeNewtonRoot(const ADReal & x1, const ADReal x2)
{
  // Machine floating-point precision
  const Real eps = std::numeric_limits<Real>::epsilon();
  // Declare some stuff here
  ADReal xh, xl;
  ADReal fl = iterativeResidual(x1);
  ADReal fh = iterativeResidual(x2);
  if (fl * fh > 0.0) // fl and fh are of the same sign
    throw MooseException("LynxCreepModel: in 'safeNewtonRoot', the points x1 and x2 "
                         "must bracket a root of the function!\nx1 = ",
                         x1,
                         " and f(x1) = ",
                         fl,
                         "\nx2 = ",
                         x2,
                         " and f(x2) = ",
                         fh,
                         "\n");

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