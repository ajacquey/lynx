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

#include "LynxDeformationBase.h"
#include "libmesh/utility.h"
#include "Assembly.h"

template <>
InputParameters
validParams<LynxDeformationBase>()
{
  InputParameters params = validParams<LynxMaterialBase>();
  params.addClassDescription("Base class calculating the deformation of a material.");
  // Coupled variables
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system.");
  params.addCoupledVar("temperature", "The temperature variable.");
  params.addCoupledVar("lithostatic_pressure", "The lithostatic pressure variable.");
  params.addCoupledVar("dynamic_pressure", "The dynamic pressure variable.");
  // Strain parameters
  params.addParam<MooseEnum>("strain_model",
                             LynxDeformationBase::strainModel() = "small",
                             "The model to use to calculate the strain rate tensor.");
  params.addParam<bool>(
      "volumetric_locking_correction", false, "Flag to correct volumetric locking");
  params.suppressParameter<bool>("use_displaced_mesh");
  // Elastic moduli parameters
  params.addRangeCheckedParam<std::vector<Real>>(
      "bulk_modulus", "bulk_modulus >= 0.0", "The drained bulk modulus of the material.");
  params.addRangeCheckedParam<std::vector<Real>>(
      "shear_modulus", "shear_modulus >= 0.0", "The shear modulus of the material.");
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
  params.addParam<MooseEnum>(
      "viscous_update",
      LynxDeformation::viscousUpdate() = "direct",
      "The type of update for the viscous strain (direct, Brent or Newton method).");
  params.addParam<Real>("tolerance", 1.0e-20, "Tolerance for the root finding algorithms.");
  params.addRangeCheckedParam<unsigned int>(
      "max_iterations",
      500,
      "max_iterations>0",
      "Maximum number of iterations in the root finding algorithms.");
  params.addParam<std::vector<Real>>(
      "initial_viscosity", "The vector of initial viscosity for purely viscous deformation.");
  params.addParam<Real>("background_strain_rate",
                        "The background strain rate for purely viscous deformation.");
  params.addParam<std::vector<Real>>(
      "eta_min", "The lower threshold for the effective viscosity for purely viscous deformation.");
  params.addParam<std::vector<Real>>(
      "eta_max", "The upper threshold for the effective viscosity for purely viscous deformation.");
  return params;
}

LynxDeformationBase::LynxDeformationBase(const InputParameters & parameters)
  : LynxMaterialBase(parameters),
    // Coupled variables
    _ndisp(coupledComponents("displacements")),
    _grad_disp(3),
    _grad_disp_old(3),
    _coupled_temp(isCoupled("temperature")),
    _temp(_coupled_temp ? coupledValue("temperature") : _zero),
    _coupled_plith(isCoupled("lithostatic_pressure")),
    _plith(_coupled_plith ? coupledValue("lithostatic_pressure") : _zero),
    _plith_old(_coupled_plith ? coupledValueOld("lithostatic_pressure") : _zero),
    _coupled_pdyn(isCoupled("dynamic_pressure")),
    _pdyn(_coupled_pdyn ? coupledValue("dynamic_pressure") : _zero),
    // Strain parameters
    _strain_model(getParam<MooseEnum>("strain_model")),
    _vol_locking_correction(getParam<bool>("volumetric_locking_correction")),
    _deformation_gradient(_fe_problem.getMaxQps()),
    _deformation_gradient_old(_fe_problem.getMaxQps()),
    _current_elem_volume(_assembly.elemVolume()),
    // Elastic moduli parameters
    _has_bulk_modulus(isParamValid("bulk_modulus")),
    _has_shear_modulus(isParamValid("shear_modulus")),
    _bulk_modulus(_has_bulk_modulus ? getLynxParam<Real>("bulk_modulus")
                                    : std::vector<Real>(_n_composition, 0.0)),
    _shear_modulus(_has_shear_modulus ? getLynxParam<Real>("shear_modulus")
                                      : std::vector<Real>(_n_composition, 0.0)),
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
                                          : std::vector<Real>(0.0)),
    _E_dislocation((_has_dislocation_creep && isParamValid("E_dislocation"))
                       ? getLynxParam<Real>("E_dislocation")
                       : std::vector<Real>(_n_composition, 0.0)),
    _V_dislocation((_has_dislocation_creep && isParamValid("V_dislocation"))
                       ? getLynxParam<Real>("V_dislocation")
                       : std::vector<Real>(_n_composition, 0.0)),
    _n_dislocation((_has_dislocation_creep && isParamValid("n_dislocation"))
                       ? getLynxParam<Real>("n_dislocation")
                       : std::vector<Real>(_n_composition, 1.0)),
    _gas_constant(getParam<Real>("gas_constant")),
    _viscous_update(getParam<MooseEnum>("viscous_update")),
    _tol(getParam<Real>("tolerance")),
    _itmax(getParam<unsigned int>("max_iterations")),
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
    // Rheology boolean
    _has_elastic(_has_shear_modulus && (_shear_modulus != std::vector<Real>(_n_composition, 0.0))),
    _has_viscous(
        (_has_diffusion_creep && (_A_diffusion != std::vector<Real>(_n_composition, 0.0))) ||
        (_has_dislocation_creep && (_A_dislocation != std::vector<Real>(_n_composition, 0.0)))),
    _has_plasticity(false),
    // Strain properties
    _elastic_strain(declareProperty<RankTwoTensor>("elastic_strain")),
    _elastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("elastic_strain")),
    _strain_increment(declareProperty<RankTwoTensor>("strain_increment")),
    _spin_tensor(declareProperty<RankTwoTensor>("spin_tensor")),
    _thermal_strain_incr(getDefaultMaterialProperty<RankTwoTensor>("thermal_strain_increment")),
    // Viscous properties
    _viscous_strain_incr(declareProperty<RankTwoTensor>("viscous_strain_increment")),
    _eta_eff(declareProperty<Real>("effective_viscosity")),
    // Plastic properties
    _plastic_strain_incr(declareProperty<RankTwoTensor>("plastic_strain_increment")),
    _plastic_yield_function(declareProperty<Real>("plastic_yield_function")),
    // Stress properties
    _stress(declareProperty<RankTwoTensor>("stress")),
    _K(declareProperty<Real>("bulk_modulus")),
    _G(declareProperty<Real>("shear_modulus")),
    _tangent_modulus(declareProperty<RankFourTensor>("tangent_modulus")),
    _inelastic_heat(declareProperty<Real>("inelastic_heat")),
    _dinelastic_heat_dstrain(declareProperty<RankTwoTensor>("dinelastic_heat_dstrain"))
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError(
        "The number of variables supplied in 'displacements' must match the mesh dimension.");

  // Fetch coupled variables and gradients (as stateful properties if necessary)
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _grad_disp[i] = &coupledGradient("displacements", i);
    if (_fe_problem.isTransient())
      _grad_disp_old[i] = &coupledGradientOld("displacements", i);
    else
      _grad_disp_old[i] = &_grad_zero;
  }

  // Set unused dimensions to zero
  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _grad_disp[i] = &_grad_zero;
    _grad_disp_old[i] = &_grad_zero;
  }

  // Creep structures
  if (_has_diffusion_creep)
    _diffusion_creep = new diffusion_creep();
  if (_has_dislocation_creep)
    _dislocation_creep = new dislocation_creep();
  if (_viscous_update != "direct")
    _iterative_viscous = new iterative_viscous();

  // Some RankFourTensor utilities
  _identity_two = RankTwoTensor(RankTwoTensor::initIdentity);
  _volumetric_four = _identity_two.outerProduct(_identity_two);
  _identity_four = RankFourTensor(RankFourTensor::initIdentitySymmetricFour);
  _deviatoric_four = _identity_four - _volumetric_four / 3.0;
}

MooseEnum
LynxDeformationBase::strainModel()
{
  return MooseEnum("small=0 finite=1");
}

MooseEnum
LynxDeformationBase::viscousUpdate()
{
  return MooseEnum("direct=1 brent=2 newton=3");
}

void
LynxDeformationBase::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
}

void
LynxDeformationBase::computeProperties()
{
  computeStrainIncrement();
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
    computeQpProperties();
}

void
LynxDeformationBase::computeQpProperties()
{
  computeQpDeformation();
  computeQpThermalSources();
}

void
LynxDeformationBase::computeStrainIncrement()
{
  Real vol_strain_incr = 0.0;
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    RankTwoTensor grad_tensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);
    RankTwoTensor grad_tensor_old(
        (*_grad_disp_old[0])[_qp], (*_grad_disp_old[1])[_qp], (*_grad_disp_old[2])[_qp]);

    _deformation_gradient[_qp] = grad_tensor;
    _deformation_gradient[_qp].addIa(1.0);

    _deformation_gradient_old[_qp] = grad_tensor_old;
    _deformation_gradient_old[_qp].addIa(1.0);

    switch (_strain_model)
    {
      case 0: // SMALL STRAIN
        calculateSmallStrain(grad_tensor, grad_tensor_old);
        break;
      case 1: // FINITE STRAIN
        calculateFiniteStrain(grad_tensor, grad_tensor_old);
        break;
      default:
        mooseError("Unknown strain model. Specify 'small' or 'finite'!");
    }

    // Thermal strain correction
    _strain_increment[_qp] -= _thermal_strain_incr[_qp];

    if (_vol_locking_correction)
      vol_strain_incr += _strain_increment[_qp].trace() * _JxW[_qp] * _coord[_qp];
  }

  if (_vol_locking_correction)
  {
    vol_strain_incr /= _current_elem_volume;

    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
    {
      Real trace = _strain_increment[_qp].trace();
      _strain_increment[_qp](0, 0) += (vol_strain_incr - trace) / 3.0;
      _strain_increment[_qp](1, 1) += (vol_strain_incr - trace) / 3.0;
      _strain_increment[_qp](2, 2) += (vol_strain_incr - trace) / 3.0;
    }
  }
}

void
LynxDeformationBase::calculateSmallStrain(const RankTwoTensor & grad_tensor,
                                          const RankTwoTensor & grad_tensor_old)
{
  RankTwoTensor A = grad_tensor - grad_tensor_old;

  _strain_increment[_qp] = 0.5 * (A + A.transpose());
  _spin_tensor[_qp] = 0.5 * (A - A.transpose());
}

void
LynxDeformationBase::calculateFiniteStrain(const RankTwoTensor & grad_tensor,
                                           const RankTwoTensor & grad_tensor_old)
{
  RankTwoTensor F = grad_tensor;
  RankTwoTensor F_old = grad_tensor_old;
  F.addIa(1.0);
  F_old.addIa(1.0);

  // Increment gradient
  RankTwoTensor L = -F_old * F.inverse();
  L.addIa(1.0);

  _strain_increment[_qp] = 0.5 * (L + L.transpose());
  _spin_tensor[_qp] = 0.5 * (L - L.transpose());
}

void
LynxDeformationBase::computeQpDeformation()
{
  // Initialize deformation
  initializeQpDeformation();

  // Update elastic moduli
  elasticModuli();

  // Update the volumetric part of the deformation
  Real pressure = volumetricDeformation();

  // Update the deviatoric part of the deformation
  RankTwoTensor stress_dev = deviatoricDeformation(pressure);

  // Plastic correction
  if (_has_plasticity && (_G[_qp] != 0.0) && (_K[_qp] != 0.0))
    plasticCorrection(pressure, stress_dev);

  // Form the total stress tensor
  reformStressTensor(pressure, stress_dev);

  // Additional correction
  if (_has_plasticity && (_G[_qp] != 0.0) && (_K[_qp] != 0.0))
    damageCorrection();

  // Update tangent operator modulus
  if (_fe_problem.currentlyComputingJacobian())
    tangentOperator();
}

void
LynxDeformationBase::initializeQpDeformation()
{
  // Initialze elastic strain
  _elastic_strain[_qp] = spinRotation(_elastic_strain_old[_qp]) + _strain_increment[_qp];

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
LynxDeformationBase::elasticModuli()
{
  // Bulk modulus
  _K[_qp] = _has_bulk_modulus ? averageProperty(_bulk_modulus) : 0.0;

  // Shear modulus
  _G[_qp] = _has_shear_modulus ? averageProperty(_shear_modulus) : 0.0;
}

Real
LynxDeformationBase::volumetricDeformation()
{
  if (_coupled_pdyn)
    return _plith[_qp] + _pdyn[_qp];
  else
  {
    Real pressure = -_K[_qp] * _elastic_strain[_qp].trace();
    // Lithostatic pressure
    pressure += _plith[_qp];
    return pressure;
  }
}

RankTwoTensor
LynxDeformationBase::deviatoricDeformation(const Real & pressure)
{
  // Update if elasticity is provided
  if (_G[_qp] != 0.0)
  {
    RankTwoTensor stress_dev = 2.0 * _G[_qp] * _elastic_strain[_qp].deviatoric();
    Real eqv_stress = std::sqrt(1.5) * stress_dev.L2norm();
    RankTwoTensor flow_dir = (eqv_stress != 0.0) ? stress_dev / eqv_stress : RankTwoTensor();

    // Viscous correction
    Real eqv_v_strain_incr = viscousIncrement(pressure, eqv_stress);

    _viscous_strain_incr[_qp] = 1.5 * eqv_v_strain_incr * flow_dir;
    stress_dev -= 3.0 * _G[_qp] * eqv_v_strain_incr * flow_dir;
    _elastic_strain[_qp] -= _viscous_strain_incr[_qp];

    return stress_dev;
  }
  else
  {
    if (!_fe_problem.isTransient())
      mooseError("Cannot run a purely viscous model in steady state!");
    if (!_has_viscous)
      mooseError("Shear modulus and initial viscosity are null! Cannot compute rheology!\n");

    _viscous_strain_incr[_qp] = _strain_increment[_qp].deviatoric();
    _eta_eff[_qp] = computeStokeEffectiveViscosity(pressure);
    // Here calculate stress_dev based on Stoke flow.
    return 2.0 * _eta_eff[_qp] * _strain_increment[_qp].deviatoric() / _dt;
  }
}

void
LynxDeformationBase::reformStressTensor(const Real & pressure, const RankTwoTensor & stress_dev)
{
  _stress[_qp] = stress_dev;
  _stress[_qp].addIa(-pressure);
}

void
LynxDeformationBase::damageCorrection()
{
}

void
LynxDeformationBase::tangentOperator()
{
  // Build flow_direction_dyad
  RankTwoTensor stress_dev = _stress[_qp].deviatoric();
  Real eqv_stress = std::sqrt(1.5) * stress_dev.L2norm();
  const RankTwoTensor flow_direction =
      (eqv_stress != 0.0) ? stress_dev / eqv_stress : RankTwoTensor();
  const RankFourTensor flow_direction_dyad = flow_direction.outerProduct(flow_direction);

  // Elastic operator
  _tangent_modulus[_qp] = (_G[_qp] != 0.0) ? 2.0 * _G[_qp] * _deviatoric_four
                                           : 2.0 * _eta_eff[_qp] * _deviatoric_four / _dt;

  _tangent_modulus[_qp] += _K[_qp] * _volumetric_four;

  // Viscous correction
  RankFourTensor viscous_operator = viscousTangentOperator(flow_direction_dyad);
  _tangent_modulus[_qp] -= 2.0 * _G[_qp] * viscous_operator;

  // Plastic correction
  RankFourTensor plastic_operator = plasticTangentOperator(flow_direction, flow_direction_dyad);
  _tangent_modulus[_qp] = _tangent_modulus[_qp] * (_identity_four - plastic_operator);

  // Additional correction
  RankFourTensor tme = (_identity_four - viscous_operator) * (_identity_four - plastic_operator);
  RankFourTensor damage_operator = damageTangentOperator(tme);
  _tangent_modulus[_qp] -= damage_operator;

  // Spin correction
  updateSpinTangentModulus();

  // Update tangent modulus for finite strain
  if (_strain_model == 1)
    finiteTangentOperator();
}

RankFourTensor
LynxDeformationBase::viscousTangentOperator(const RankFourTensor & flow_direction_dyad)
{
  if ((_eta_eff[_qp] != 0.0) && (_G[_qp] != 0.0))
    return (1.0 - _stress_corr_v) * _deviatoric_four +
           1.5 * (_stress_corr_v - _dq_dq_tr_v) * flow_direction_dyad;
  else
    return RankFourTensor();
}

RankFourTensor
LynxDeformationBase::plasticTangentOperator(const RankTwoTensor & flow_direction,
                                            const RankFourTensor & flow_direction_dyad)
{
  if (_has_plasticity && (_G[_qp] != 0.0) && (_K[_qp] != 0.0))
  {
    RankFourTensor plastic_operator = (1.0 - _stress_corr_p) * _deviatoric_four +
                                      1.5 * (_stress_corr_p - _dq_dq_tr_p) * flow_direction_dyad;
    plastic_operator +=
        0.5 * _K[_qp] / _G[_qp] * _dq_dp_tr_p * flow_direction.outerProduct(_identity_two);
    plastic_operator += (1.0 - _dp_dp_tr_p) * _volumetric_four / 3.0;
    plastic_operator +=
        3.0 * _G[_qp] / _K[_qp] * _dp_dq_tr_p * _identity_two.outerProduct(flow_direction) / 3.0;

    return plastic_operator;
  }
  else
    return RankFourTensor();
}

RankFourTensor
LynxDeformationBase::damageTangentOperator(const RankFourTensor & /*tme*/)
{
  return RankFourTensor();
}

void
LynxDeformationBase::finiteTangentOperator()
{
  RankTwoTensor Id = RankTwoTensor(RankTwoTensor::initIdentity);
  RankFourTensor dL_dFinv = -_deformation_gradient_old[_qp].mixedProductIkJl(Id);
  RankTwoTensor dfr_gdt_inv = _deformation_gradient[_qp].inverse();
  RankFourTensor dFinv_dF = -dfr_gdt_inv.mixedProductIkJl(dfr_gdt_inv.transpose());

  RankFourTensor dL_dgrad_u = dL_dFinv * dFinv_dF;

  _tangent_modulus[_qp] = _tangent_modulus[_qp] * dL_dgrad_u;
}

Real
LynxDeformationBase::computeStokeEffectiveViscosity(const Real & pressure)
{
  computeCreepProperties(pressure);

  Real strain_rate_II = std::sqrt(0.5) * _strain_increment[_qp].deviatoric().L2norm() / _dt;

  if (_t_step <= 1 && strain_rate_II == 0.0 && _has_initial_viscosity)
    return std::min(std::max(averageProperty(_initial_viscosity), averageProperty(_eta_min)),
                    averageProperty(_eta_max));

  strain_rate_II = (_has_background_strain_rate && _t_step <= 1 && strain_rate_II == 0.0)
                       ? _background_strain_rate
                       : strain_rate_II;

  Real one_on_eta_diff = _has_diffusion_creep ? _diffusion_creep->oneOnViscosity() : 0.0;
  Real one_on_eta_disl =
      _has_dislocation_creep ? _dislocation_creep->oneOnViscosity(strain_rate_II) : 0.0;

  Real eta = 1.0 / (one_on_eta_diff + one_on_eta_disl);

  return std::min(std::max(eta, averageProperty(_eta_min)), averageProperty(_eta_max));
}

RankTwoTensor
LynxDeformationBase::computeStokeEffectiveViscosityDerivative()
{
  RankTwoTensor strain_rate = _strain_increment[_qp].deviatoric() / _dt;
  Real strain_rate_II = std::sqrt(0.5) * strain_rate.L2norm();

  Real deta_dstrain_rate_II = _dislocation_creep->oneOnViscosityDerivative(strain_rate_II);
  RankTwoTensor dstrain_rate_II =
      (strain_rate_II != 0.0) ? 0.5 * strain_rate / strain_rate_II : RankTwoTensor();

  return -Utility::pow<2>(_eta_eff[_qp]) * deta_dstrain_rate_II * dstrain_rate_II;
}

Real
LynxDeformationBase::computeEffectiveViscosity(const Real & pressure, const Real & eqv_stress)
{
  computeCreepProperties(pressure);

  Real eff_creep_rate = 0.0;
  Real eta_eff = 0.0;

  if (eqv_stress != 0.0)
  {
    if (_viscous_update == "direct")
    {
      Real diff_creep_rate = 0.0, disl_creep_rate = 0.0;
      if (_has_diffusion_creep)
        diff_creep_rate = _diffusion_creep->creepRate(eqv_stress);
      if (_has_dislocation_creep)
        disl_creep_rate = _dislocation_creep->creepRate(eqv_stress);

      eff_creep_rate = diff_creep_rate + disl_creep_rate;
      eta_eff = (eff_creep_rate != 0.0) ? eqv_stress / (3.0 * eff_creep_rate) : 0.0;
    }
    else // Iterative method
    {
      _iterative_viscous->fill(eqv_stress, _G[_qp], _dt);
      if (_has_diffusion_creep)
        _iterative_viscous->fillDiff(_diffusion_creep);
      if (_has_dislocation_creep)
        _iterative_viscous->fillDisl(_dislocation_creep);

      // Iterative methods require two starting points bracketing the root
      Real x_pos = 0.0;
      Real x_neg = eqv_stress / (3.0 * _G[_qp] * _dt);

      if (_iterative_viscous->func(x_pos) == 0.0)
        eff_creep_rate = x_pos;
      else if (_iterative_viscous->func(x_neg) == 0.0)
        eff_creep_rate = x_neg;
      else if (_iterative_viscous->func(x_pos) * _iterative_viscous->func(x_neg) < 0.0)
      {
        switch (_viscous_update)
        {
          case 2:
            eff_creep_rate = rootBrent(*_iterative_viscous, x_pos, x_neg);
            break;
          case 3:
            eff_creep_rate = rootNewtonSafe(*_iterative_viscous, x_pos, x_neg);
            break;
          default:
            mooseError("LynxDeformation: unknown update procedure type!");
        }
      }
      else
        throw MooseException("LynxDeformation: could not fint bracketing values of creep strain "
                             "rate for iterative methods!\nx_pos is " +
                             Moose::stringify(x_pos) + " and F(x_pos) is " +
                             Moose::stringify(_iterative_viscous->func(x_pos)) + "\nx_neg is " +
                             Moose::stringify(x_neg) + " and F(x_neg) is " +
                             Moose::stringify(_iterative_viscous->func(x_neg)) + "\n");

      eta_eff = (eff_creep_rate != 0.0)
                    ? (eqv_stress - 3.0 * _G[_qp] * eff_creep_rate * _dt) / (3.0 * eff_creep_rate)
                    : 0.0;
    }
  }

  return std::min(std::max(eta_eff, averageProperty(_eta_min)), averageProperty(_eta_max));
}

Real
LynxDeformationBase::computeCreepRate(const Real A, const Real n, const Real eqv_stress)
{
  Real s_II = eqv_stress / std::sqrt(3.0);
  return 2.0 / std::sqrt(3.0) * A * std::pow(s_II, n);
}

Real
LynxDeformationBase::viscousIncrement(const Real & pressure, const Real & eqv_stress)
{
  // Compute effective viscosity
  _eta_eff[_qp] = computeEffectiveViscosity(pressure, eqv_stress);

  Real eqv_v_strain_incr = 0.0;

  if (_eta_eff[_qp] != 0.0) // viscous correction
  {
    Real maxwell_time = _eta_eff[_qp] / _G[_qp];
    eqv_v_strain_incr =
        (_dt / maxwell_time) / (1.0 + (_dt / maxwell_time)) * eqv_stress / (3.0 * _G[_qp]);

    // Tangent Operator
    if (_fe_problem.currentlyComputingJacobian())
    {
      _stress_corr_v =
          (eqv_stress != 0.0) ? (eqv_stress - 3.0 * _G[_qp] * eqv_v_strain_incr) / eqv_stress : 1.0;
      _dq_dq_tr_v = 1.0 / (1.0 + _dt / maxwell_time);
    }
  }

  return eqv_v_strain_incr;
}

RankTwoTensor
LynxDeformationBase::spinRotation(const RankTwoTensor & tensor)
{
  return tensor + _spin_tensor[_qp] * tensor.deviatoric() - tensor.deviatoric() * _spin_tensor[_qp];
}

void
LynxDeformationBase::updateSpinTangentModulus()
{
  RankTwoTensor strain_el_dev_old = _elastic_strain_old[_qp].deviatoric();
  RankTwoTensor Id = RankTwoTensor(RankTwoTensor::initIdentity);

  if (_G[_qp] != 0.0)
    _tangent_modulus[_qp] +=
        _G[_qp] * (Id.mixedProductIkJl(strain_el_dev_old.transpose()) -
                   Id.mixedProductIlJk(strain_el_dev_old.transpose())) -
        _G[_qp] * (strain_el_dev_old.mixedProductIkJl(Id) - strain_el_dev_old.mixedProductIlJk(Id));
}

void
LynxDeformationBase::computeCreepProperties(const Real & pressure)
{
  updateCreepParameters();
  if (_coupled_temp)
  {
    Real Q_diff = 0.0, Q_disl = 0.0;
    Real RT = _gas_constant * _temp[_qp];

    if (_has_diffusion_creep)
    {
      Q_diff = (_diffusion_creep->_act_energy + _diffusion_creep->_act_vol * pressure) / RT;
      _diffusion_creep->_A = _diffusion_creep->_A0 * std::exp(-Q_diff);
    }

    if (_has_dislocation_creep)
    {
      Q_disl = (_dislocation_creep->_act_energy + _dislocation_creep->_act_vol * pressure) / RT;
      _dislocation_creep->_A = _dislocation_creep->_A0 * std::exp(-Q_disl);
    }
  }
  else
  {
    if (_has_diffusion_creep)
      _diffusion_creep->_A = _diffusion_creep->_A0;
    if (_has_dislocation_creep)
      _dislocation_creep->_A = _dislocation_creep->_A0;
  }
}

void
LynxDeformationBase::updateCreepParameters()
{
  if (_has_diffusion_creep)
    _diffusion_creep->fill(averageProperty(_A_diffusion),
                           averageProperty(_E_diffusion),
                           averageProperty(_V_diffusion));

  if (_has_dislocation_creep)
    _dislocation_creep->fill(averageProperty(_A_dislocation),
                             averageProperty(_n_dislocation),
                             averageProperty(_E_dislocation),
                             averageProperty(_V_dislocation));
}

// Finding root of function using the Van Wijngaarden-Dekker-Brent method or Brent method in short.
// The values x1 and x2 must bracket the root (f(x1) * f(x2) < 0)
// See 9.3 of Numerical Recipes. The Art of Scientific Computing, 3rd Edition, 2007, (C++ code)
Real
LynxDeformationBase::rootBrent(iterative_viscous & viscous_model, const Real x1, const Real x2)
{
  // Machine floating-point precision
  const Real eps = std::numeric_limits<Real>::epsilon();
  // Declare some stuff here
  Real a = x1, b = x2, c = x2, d, e, fa = viscous_model.func(a), fb = viscous_model.func(b), fc, p,
       q, r, s, tol1, xm;
  // Chek if x1 and x2 bracket a root
  if (fa * fb > 0.0) // fa and fb are of the same sign
    throw MooseException("LynxViscousIterativeUtils: in Brent's method, the points x1 and x2 must "
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
    fb = viscous_model.func(b);
  }
  throw MooseException(
      "LynxViscousIterativeUtils: maximum number of iterations exceeded in solveBrent!");
}

// Finding root of function using a safe Newton-Raphson method (combinsaison of Newton-Raphson and
// bisection). The values x1 and x2 must bracket the root (f(x1) * f(x2) < 0) See 9.4 of Numerical
// Recipes. The Art of Scientific Computing, 3rd Edition, 2007, (C++ code)
Real
LynxDeformationBase::rootNewtonSafe(iterative_viscous & viscous_model, const Real x1, const Real x2)
{
  // Machine floating-point precision
  const Real eps = std::numeric_limits<Real>::epsilon();
  // Declare some stuff here
  Real xh, xl;
  Real fl = viscous_model.func(x1);
  Real fh = viscous_model.func(x2);
  if (fl * fh > 0.0) // fl and fh are of the same sign
    throw MooseException("LynxDeformation: in safe Newton-Raphson's method, the points x1 and x2 "
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
  Real rts = 0.5 * (x1 + x2);
  // The "stepsize before last"
  Real dxold = std::abs(x2 - x1);
  // and the last step
  Real dx = dxold;
  Real f = viscous_model.func(rts);
  Real df = viscous_model.deriv(rts);
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
      Real temp = rts;
      rts -= dx;
      // Change in root is negligible.
      if (temp == rts)
        return rts;
    }
    // Convergence criterion
    if (std::abs(dx) < (2.0 * eps * std::abs(rts) + 0.5 * _tol))
      return rts;
    f = viscous_model.func(rts);
    df = viscous_model.deriv(rts);
    // The one new function evaluation per iteration.
    // Maintain the bracket on the root
    if (f < 0.0)
      xl = rts;
    else
      xh = rts;
  }
  throw MooseException(
      "LynxDeformation: maximum number of iterations exceeded in solveNewtonSafe!");
}

void
LynxDeformationBase::computeQpThermalSources()
{
  RankTwoTensor inelastic_strain_incr = _viscous_strain_incr[_qp] + _plastic_strain_incr[_qp] + _thermal_strain_incr[_qp];
  _inelastic_heat[_qp] = _stress[_qp].doubleContraction(inelastic_strain_incr) / _dt;

  // if (_fe_problem.currentlyComputingJacobian())
  // {
  //   if (inelastic_strain_incr == RankTwoTensor())
  //     // no inelastic deformation, so _elasticity_tensor = _tangent_modulus
  //     _dinelastic_heat_dstrain[_qp] = RankTwoTensor();
  //   else
  //   {
  //     // Build deviatoric tensor for tangent modulus
  //     RankTwoTensor identity_two(RankTwoTensor::initIdentity);
  //     RankFourTensor deviatoric_four(RankFourTensor::initIdentitySymmetricFour);
  //     deviatoric_four -= identity_two.outerProduct(identity_two) / 3.0;

  //     _dinelastic_heat_dstrain[_qp] =
  //         inelastic_strain_incr.initialContraction(_tangent_modulus[_qp]);
  //     _dinelastic_heat_dstrain[_qp] += _stress[_qp];
  //     _dinelastic_heat_dstrain[_qp] -= _stress[_qp].deviatoric().initialContraction(
  //         _tangent_modulus[_qp] * deviatoric_four / (2.0 * _G[_qp]));
  //     _dinelastic_heat_dstrain[_qp].addIa(-_stress[_qp].trace() * _tangent_modulus[_qp].sum3x3()
  //     /
  //                                         (3.0 * _K[_qp]));
  //     _dinelastic_heat_dstrain[_qp] /= _dt;
  //   }
  // }
}
