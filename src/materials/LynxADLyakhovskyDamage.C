// /******************************************************************************/
// /*                            This file is part of                            */
// /*                       LYNX, a MOOSE-based application                      */
// /*                    Lithosphere dYnamic Numerical toolboX                   */
// /*                                                                            */
// /*          Copyright (C) 2017 by Antoine B. Jacquey and Mauro Cacace         */
// /*             GFZ Potsdam, German Research Centre for Geosciences            */
// /*                                                                            */
// /*                Licensed under GNU General Public License 3,                */
// /*                       please see LICENSE for details                       */
// /*                  or http://www.gnu.org/licenses/gpl.html                   */
// /******************************************************************************/

// #include "LynxADLyakhovskyDamage.h"
// #include "ElasticityTensorTools.h"

// registerMooseObject("LynxApp", LynxADLyakhovskyDamage);

// InputParameters
// LynxADLyakhovskyDamage::validParams()
// {
//   InputParameters params = LynxADDamageModelBase::validParams();
//   params.addClassDescription("Damage rheology of Lyakhovsky.");
//   // Damage properties
//   params.addRequiredParam<std::vector<Real>>("damage_modulus",
//                                              "The third elastic modulus for damage rheology.");
//   // params.addRequiredParam<std::vector<Real>>("critical_strain_ratio", "The critical strain ratio
//   // for damage rheology.");
//   params.addRequiredParam<std::vector<Real>>(
//       "friction_angle",
//       "The friction angle of the material for the pressure-dependant part of the yield stress.");
//   params.addParam<std::vector<Real>>("cohesion",
//                                      "The constant coefficient of the yield stress corresponding "
//                                      "to the cohesion of the material.");
//   params.addRequiredParam<std::vector<Real>>("plastic_viscosity", "The plastic viscosity.");
//   return params;
// }

// LynxADLyakhovskyDamage::LynxADLyakhovskyDamage(const InputParameters & parameters)
//   : LynxADDamageModelBase(parameters),
//     _damage_modulus(getLynxParam<Real>("damage_modulus")),
//     _friction_angle(getLynxParam<Real>("friction_angle")),
//     _cohesion(isParamValid("cohesion") ? getLynxParam<Real>("cohesion")
//                                        : std::vector<Real>(_n_composition, 0.0)),
//     _plastic_viscosity(getLynxParam<Real>("plastic_viscosity"))
// {
// }

// void
// LynxADLyakhovskyDamage::elasticGuess(ADRankTwoTensor & stress,
//                                      const ADRankTwoTensor & stress_old,
//                                      const ADRankFourTensor & Cijkl,
//                                      const ADRankTwoTensor & elastic_strain_old,
//                                      const ADRankTwoTensor & elastic_strain_incr)
// {
//   // Compute Stuff for update
//   _Gam0 = averageProperty(_damage_modulus);
//   _eta_p = averageProperty(_plastic_viscosity);
//   _K = ElasticityTensorTools::getIsotropicBulkModulus(Cijkl);
//   _G = ElasticityTensorTools::getIsotropicShearModulus(Cijkl);
//   Real sin_phi = std::sin(averageProperty(_friction_angle) * libMesh::pi / 180.0);
//   _xi0 = -std::sqrt(3.0) / std::sqrt(1.0 + 1.5 * Utility::pow<2>(_K / _G * sin_phi));
//   _k = std::sqrt(3.0) * averageProperty(_cohesion) *
//        std::cos(averageProperty(_friction_angle) * libMesh::pi / 180.0);
//   _e_norm = elastic_strain_old.L2norm();
//   _xi = strainRatio(elastic_strain_old);
//   _e = (_e_norm != 0.0) ? elastic_strain_old / _e_norm : ADRankTwoTensor();
//   _mu0 = 3.0 * _G * std::sqrt(2.0 / 3.0 * (1.0 - Utility::pow<2>(_xi0) / 3.0)) / (-_K * _xi0);

//   // Stress
//   _Dijkl = damageStiffness(Cijkl);
//   stress = stress_old + _Dijkl * elastic_strain_incr;

//   // Damage force
//   ADRankTwoTensor elastic_strain_tr = elastic_strain_old + elastic_strain_incr;
//   ADReal e_norm_tr = elastic_strain_tr.L2norm();
//   ADReal xi_tr = strainRatio(elastic_strain_tr);
//   _damage_drive[_qp] = _Gam0 * Utility::pow<2>(e_norm_tr) * (xi_tr - _xi0);
// }

// ADRankFourTensor
// LynxADLyakhovskyDamage::damageStiffness(const ADRankFourTensor & Cijkl)
// {
//   ADRankTwoTensor I = ADRankTwoTensor(ADRankTwoTensor::initIdentity);
//   ADRankFourTensor deviatoric_four = ADRankFourTensor(ADRankFourTensor::initIdentitySymmetricFour);
//   -I.outerProduct(I) / 3.0;

//   // Damage correction
//   ADRankFourTensor D = Cijkl;
//   D -= _damage_old[_qp] * _Gam0 * (_xi - 2.0 * _xi0) * deviatoric_four;
//   D -= _damage_old[_qp] * _Gam0 *
//        (_e.outerProduct(I) + I.outerProduct(_e) - _xi * _e.outerProduct(_e));

//   return D;
// }

// void
// LynxADLyakhovskyDamage::damageUpdate(ADRankTwoTensor & stress,
//                                      ADRankTwoTensor & elastic_strain_incr)
// {
//   // Trial state
//   _stress_tr = stress;
//   _pressure_tr = -_stress_tr.trace() / 3.0;
//   _dev_stress_tr = std::sqrt(1.5) * _stress_tr.deviatoric().L2norm();
//   _damage_force_tr = _damage_drive[_qp];

//   // Initialize plastic strain increment
//   _plastic_strain_incr[_qp].zero();
//   _damage_incr[_qp] = 0.0;
//   _damage[_qp] = _damage_old[_qp];
//   _damage_poro_mech[_qp] = computeDamagePoroMech(0.0, 0.0, elastic_strain_incr);

//   // Return if 0 dev stress
//   if (MooseUtils::absoluteFuzzyEqual(_dev_stress_tr, 0.0))
//   {
//     _yield_function[_qp] = -1.0;
//     return;
//   }

//   // Check yield function
//   // Update yield parameters
//   updateYieldParameters(0.0, 0.0);
//   _yield_function[_qp] = yieldFunction();
//   if (_yield_function[_qp] <= _abs_tol) // Elastic
//     return;

//   // Viscoplastic update
//   ADReal gamma_s = 0.0, gamma_d = 0.0;
//   returnMap(gamma_s, gamma_d);

//   // Update quantities
//   // Update yield parameters
//   updateYieldParameters(gamma_s, gamma_d);
//   _yield_function[_qp] = yieldFunction();
//   _plastic_strain_incr[_qp] = reformPlasticStrainTensor(gamma_s);
//   _damage_poro_mech[_qp] = computeDamagePoroMech(gamma_s, gamma_d, elastic_strain_incr);
//   elastic_strain_incr -= _plastic_strain_incr[_qp];
//   _damage_incr[_qp] = gamma_d * _dt;
//   _damage[_qp] = _damage_old[_qp] + _damage_incr[_qp];
//   _damage_drive[_qp] = _damage_force;
//   stress -= _Dijkl * _plastic_strain_incr[_qp];
//   stress -= _Gam0 * _e_norm * (_xi - 2.0 * _xi0) * _e * _damage_incr[_qp];
//   stress.addIa(-_Gam0 * _e_norm * _damage_incr[_qp]);
// }

// void
// LynxADLyakhovskyDamage::updateElasticModuli(ADReal & K, ADReal & G, const ADRankTwoTensor & elastic_strain)
// {
//   ADReal xi = strainRatio(elastic_strain);
//   ADReal e_norm = elastic_strain.L2norm();

//   K -= _damage[_qp] * _Gam0 * (2.0 / 3.0 * (xi - _xi0) + xi / 3.0 * (1.0 - Utility::pow<2>(xi) / 3.0));
//   G -= 0.5 * _damage[_qp] * _Gam0 * ((xi - 2.0 * _xi0) - xi * (1.0 - Utility::pow<2>(xi) / 3.0));
// }

// void
// LynxADLyakhovskyDamage::returnMap(ADReal & gamma_s, ADReal & gamma_d)
// {
//   // Initial residual
//   ADReal ress_ini = 0.0, resd_ini = 0.0;
//   residual(0.0, 0.0, ress_ini, resd_ini);
//   ADReal res_ini = std::sqrt(Utility::pow<2>(ress_ini) + Utility::pow<2>(resd_ini));
//   ADReal ress = ress_ini, resd = resd_ini;
//   ADReal res = res_ini;

//   // Initial jacobian
//   ADReal jacss = 0.0, jacdd = 0.0, jacsd = 0.0, jacds = 0.0;
//   jacobian(0.0, 0.0, jacss, jacdd, jacsd, jacds);

//   // Useful stuff
//   ADReal jac_full = jacss * jacdd - jacsd * jacds;
//   ADReal ress_full = jacdd * ress - jacsd * resd;
//   ADReal resd_full = jacss * resd - jacds * ress;

//   // Newton loop
//   for (unsigned int iter = 0; iter < _max_its; ++iter)
//   {
//     gamma_s -= ress_full / jac_full;
//     gamma_d -= resd_full / jac_full;

//     residual(gamma_s, gamma_d, ress, resd);
//     jacobian(gamma_s, gamma_d, jacss, jacdd, jacsd, jacds);
//     jac_full = jacss * jacdd - jacsd * jacds;
//     ress_full = jacdd * ress - jacsd * resd;
//     resd_full = jacss * resd - jacds * ress;
//     res = std::sqrt(Utility::pow<2>(ress) + Utility::pow<2>(resd));

//     // Convergence check
//     if ((std::abs(res) <= _abs_tol) || (std::abs(res / res_ini) <= _rel_tol))
//       return;
//   }
//   // mooseError("Failed returnMap!.\n");
//   throw MooseException("LynxADLyakhovskyDamage: maximum number of iterations exceeded in "
//                        "'returnMap'!\nInitial residual: ",
//                        res_ini,
//                        "\nResidual: ",
//                        res,
//                        "\n Dev strain rate: ",
//                        gamma_s,
//                        "\n Damage Rate: ",
//                        gamma_d,
//                        "\n");
// }

// void
// LynxADLyakhovskyDamage::residual(const ADReal & gamma_s,
//                                  const ADReal & gamma_d,
//                                  ADReal & ress,
//                                  ADReal & resd)
// {
//   // Update yield parameters
//   updateYieldParameters(gamma_s, gamma_d);

//   overStress(gamma_s, gamma_d, ress, resd);
//   ress -= _eta_p * gamma_s;
//   resd -= _eta_p * gamma_d;
// }

// void
// LynxADLyakhovskyDamage::jacobian(const ADReal & gamma_s,
//                                  const ADReal & gamma_d,
//                                  ADReal & jacss,
//                                  ADReal & jacdd,
//                                  ADReal & jacsd,
//                                  ADReal & jacds)
// {
//   // Update yield parameters and their derivatives
//   updateYieldParameters(gamma_s, gamma_d);
//   updateYieldParametersDerivS(gamma_s, gamma_d);
//   updateYieldParametersDerivD(gamma_s, gamma_d);

//   overStressDerivS(gamma_s, gamma_d, jacss, jacds);
//   overStressDerivD(gamma_s, gamma_d, jacsd, jacdd);
//   jacss -= _eta_p;
//   jacdd -= _eta_p;
// }

// void
// LynxADLyakhovskyDamage::overStress(const ADReal & gamma_s,
//                                    const ADReal & gamma_d,
//                                    ADReal & over_s,
//                                    ADReal & over_d)
// {
//   // Flow rate and directions
//   ADReal phi = flowRate();
//   ADReal rg = flowDirectionS();
//   ADReal ra = flowDirectionD();

//   over_s = phi * rg;
//   over_d = phi * ra;
// }

// void
// LynxADLyakhovskyDamage::overStressDerivS(const ADReal & gamma_s,
//                                          const ADReal & gamma_d,
//                                          ADReal & over_s_s,
//                                          ADReal & over_d_s)
// {
//   // Flow rate and directions and their derivatives
//   ADReal phi = flowRate();
//   ADReal rg = flowDirectionS();
//   ADReal ra = flowDirectionD();
//   ADReal dphi_dgammaS = dFlowRatedS();
//   ADReal dphi_dgammaD = dFlowRatedD();
//   ADReal drg_dgammaS = dFlowDirectionSdS();
//   ADReal drg_dgammaD = dFlowDirectionSdD();
//   ADReal dra_dgammaS = dFlowDirectionSdS();
//   ADReal dra_dgammaD = dFlowDirectionSdD();

//   over_s_s = dphi_dgammaS * rg + phi * drg_dgammaS;
//   over_d_s = dphi_dgammaS * ra + phi * dra_dgammaS;
// }

// void
// LynxADLyakhovskyDamage::overStressDerivD(const ADReal & gamma_s,
//                                          const ADReal & gamma_d,
//                                          ADReal & over_s_d,
//                                          ADReal & over_d_d)
// {
//   // Flow rate and directions and their derivatives
//   ADReal phi = flowRate();
//   ADReal rg = flowDirectionS();
//   ADReal ra = flowDirectionD();
//   ADReal dphi_dgammaS = dFlowRatedS();
//   ADReal dphi_dgammaD = dFlowRatedD();
//   ADReal drg_dgammaS = dFlowDirectionSdS();
//   ADReal drg_dgammaD = dFlowDirectionSdD();
//   ADReal dra_dgammaS = dFlowDirectionSdS();
//   ADReal dra_dgammaD = dFlowDirectionSdD();

//   over_s_d = dphi_dgammaD * rg + phi * drg_dgammaD;
//   over_d_d = dphi_dgammaD * ra + phi * dra_dgammaD;
// }

// void
// LynxADLyakhovskyDamage::updateYieldParameters(const ADReal & gamma_s, const ADReal & gamma_d)
// {
//   _pressure = _pressure_tr -
//               std::sqrt(1.5) * _damage_old[_qp] * _Gam0 *
//                   std::pow(1.0 - Utility::pow<2>(_xi) / 3.0, 1.5) * gamma_s * _dt +
//               _Gam0 * _e_norm * (1.0 + _xi / 3.0 * (_xi - 2.0 * _xi0)) * gamma_d * _dt;
//   _dev_stress = _dev_stress_tr -
//                 (3.0 * _G - 1.5 * _damage_old[_qp] * _Gam0 *
//                                 ((_xi - 2.0 * _xi0) - _xi * (1.0 - Utility::pow<2>(_xi) / 3.0))) *
//                     gamma_s * _dt -
//                 _Gam0 * (_xi - 2.0 * _xi0) * _e_norm *
//                     std::sqrt(1.5 * (1.0 - Utility::pow<2>(_xi) / 3.0)) * gamma_d * _dt;
//   _mua = (1.0 + 0.5 * (_damage_old[_qp] + gamma_d * _dt) * _Gam0 / _G * _xi0) /
//          (1.0 - (_damage_old[_qp] + gamma_d * _dt) * _Gam0 / _K *
//                     (1.0 - Utility::pow<2>(_xi0) / 3.0) / _xi0);
//   _damage_force = _damage_force_tr - _Gam0 * _e_norm * (_xi - 2.0 * _xi0) *
//                                          std::sqrt(1.5 * (1.0 - Utility::pow<2>(_xi) / 3.0)) *
//                                          gamma_s * _dt;
// }

// void
// LynxADLyakhovskyDamage::updateYieldParametersDerivS(const ADReal & gamma_s, const ADReal & gamma_d)
// {
//   _dpressure_dS = -std::sqrt(1.5) * _damage_old[_qp] * _Gam0 *
//                   std::pow(1.0 - Utility::pow<2>(_xi) / 3.0, 1.5) * _dt;
//   _ddev_stress_dS =
//       -(3.0 * _G - 1.5 * _damage_old[_qp] * _Gam0 *
//                        ((_xi - 2.0 * _xi0) - _xi * (1.0 - Utility::pow<2>(_xi) / 3.0))) *
//       _dt;
//   _dmua_dS = 0.0;
//   _ddamage_force_dS = -_Gam0 * _e_norm * (_xi - 2.0 * _xi0) *
//                       std::sqrt(1.5 * (1.0 - Utility::pow<2>(_xi) / 3.0)) * _dt;
// }

// void
// LynxADLyakhovskyDamage::updateYieldParametersDerivD(const ADReal & gamma_s, const ADReal & gamma_d)
// {
//   _dpressure_dD = _Gam0 * _e_norm * (1.0 + _xi / 3.0 * (_xi - 2.0 * _xi0)) * _dt;
//   _ddev_stress_dD = -_Gam0 * (_xi - 2.0 * _xi0) * _e_norm *
//                     std::sqrt(1.5 * (1.0 - Utility::pow<2>(_xi) / 3.0)) * _dt;
//   _dmua_dD =
//       (0.5 * _Gam0 / _G * _xi0 + _mua * _Gam0 / _K * (1.0 - Utility::pow<2>(_xi0) / 3.0) / _xi0) *
//       _dt /
//       (1.0 - (_damage_old[_qp] + gamma_d * _dt) * _Gam0 / _K * (1.0 - Utility::pow<2>(_xi0) / 3.0) /
//                  _xi0);
//   ;
//   _ddamage_force_dD = 0.0;
// }

// ADReal
// LynxADLyakhovskyDamage::flowRate()
// {
//   return yieldFunction();
// }

// ADReal
// LynxADLyakhovskyDamage::flowDirectionS()
// {
//   return dYieldFunctiondDev();
// }

// ADReal
// LynxADLyakhovskyDamage::flowDirectionD()
// {
//   return dYieldFunctiondDam();
// }

// ADReal
// LynxADLyakhovskyDamage::dFlowRatedS()
// {
//   ADReal df_ddev_stress = dYieldFunctiondDev();
//   ADReal df_dpressure = dYieldFunctiondPres();
//   ADReal df_dmua = dYieldFunctiondMua();

//   return df_ddev_stress * _ddev_stress_dS + df_dpressure * _dpressure_dS + df_dmua * _dmua_dS;
// }

// ADReal
// LynxADLyakhovskyDamage::dFlowRatedD()
// {
//   ADReal df_ddev_stress = dYieldFunctiondDev();
//   ADReal df_dpressure = dYieldFunctiondPres();
//   ADReal df_dmua = dYieldFunctiondMua();

//   return df_ddev_stress * _ddev_stress_dD + df_dpressure * _dpressure_dD + df_dmua * _dmua_dD;
// }

// ADReal
// LynxADLyakhovskyDamage::dFlowDirectionSdS()
// {
//   ADReal d2f_ddev_stress2 = d2YieldFunctiondDev2();
//   ADReal d2f_ddev_stress_dpressure = d2YieldFunctiondDevdPres();
//   ADReal d2f_ddev_stress_dmua = dYieldFunctiondMua();

//   return d2f_ddev_stress2 * _ddev_stress_dS + d2f_ddev_stress_dpressure * _dpressure_dS +
//          d2f_ddev_stress_dmua * _dmua_dS;
// }

// ADReal
// LynxADLyakhovskyDamage::dFlowDirectionSdD()
// {
//   ADReal d2f_ddev_stress2 = d2YieldFunctiondDev2();
//   ADReal d2f_ddev_stress_dpressure = d2YieldFunctiondDevdPres();
//   ADReal d2f_ddev_stress_dmua = d2YieldFunctiondDevdMua();

//   return d2f_ddev_stress2 * _ddev_stress_dD + d2f_ddev_stress_dpressure * _dpressure_dD +
//          d2f_ddev_stress_dmua * _dmua_dD;
// }

// ADReal
// LynxADLyakhovskyDamage::dFlowDirectionDdS()
// {
//   ADReal d2f_ddam_ddev_stress = d2YieldFunctiondDamdDev();
//   ADReal d2f_ddam_dpressure = d2YieldFunctiondDevdPres();
//   ADReal d2f_ddam_dmua = d2YieldFunctiondDamdMua();
//   ADReal d2f_ddam2 = d2YieldFunctiondDam2();

//   return d2f_ddam_ddev_stress * _ddev_stress_dS + d2f_ddam_dpressure * _dpressure_dS +
//          d2f_ddam_dmua * _dmua_dS + d2f_ddam2 * _ddamage_force_dS;
// }

// ADReal
// LynxADLyakhovskyDamage::dFlowDirectionDdD()
// {
//   ADReal d2f_ddam_ddev_stress = d2YieldFunctiondDamdDev();
//   ADReal d2f_ddam_dpressure = d2YieldFunctiondDevdPres();
//   ADReal d2f_ddam_dmua = d2YieldFunctiondDamdMua();
//   ADReal d2f_ddam2 = d2YieldFunctiondDam2();

//   return d2f_ddam_ddev_stress * _ddev_stress_dS + d2f_ddam_dpressure * _dpressure_dS +
//          d2f_ddam_dmua * _dmua_dS + d2f_ddam2 * _ddamage_force_dD;
// }

// ADReal
// LynxADLyakhovskyDamage::yieldFunction()
// {
//   if (MooseUtils::absoluteFuzzyEqual(_pressure, 0.0))
//     return 0.0;
//   else
//   {
//     ADReal f_squared =
//         Utility::pow<2>(_dev_stress / (_mu0 * _pressure)) + 1.0 - Utility::pow<2>(_mua) - 1.0;
//     if (MooseUtils::absoluteFuzzyEqual(f_squared, 0.0))
//       return 0.0;
//     else if (MooseUtils::absoluteFuzzyEqual(f_squared, -1.0))
//       return -1.0;
//     else
//       return std::sqrt(f_squared + 1.0) - 1.0;
//   }
// }

// ADReal
// LynxADLyakhovskyDamage::dYieldFunctiondDev()
// {
//   ADReal f = yieldFunction();

//   if (MooseUtils::absoluteFuzzyEqual(_pressure, 0.0) || MooseUtils::absoluteFuzzyEqual(f, -1.0))
//     return 0.0;
//   else
//     return _dev_stress / (Utility::pow<2>(_mu0 * _pressure) * (1.0 + f));
// }

// ADReal
// LynxADLyakhovskyDamage::dYieldFunctiondPres()
// {
//   ADReal f = yieldFunction();

//   if (MooseUtils::absoluteFuzzyEqual(_pressure, 0.0) || MooseUtils::absoluteFuzzyEqual(f, -1.0))
//     return 0.0;
//   else
//     return -Utility::pow<2>(_dev_stress) /
//            (_pressure * Utility::pow<2>(_mu0 * _pressure) * (1.0 + f));
// }

// ADReal
// LynxADLyakhovskyDamage::dYieldFunctiondDam()
// {
//   ADReal f = yieldFunction();

//   if (MooseUtils::absoluteFuzzyEqual(_damage_force, 0.0) || MooseUtils::absoluteFuzzyEqual(f, -1.0))
//     return 0.0;
//   else
//     return (1.0 - Utility::pow<2>(_mua)) / (_damage_force * (1.0 + f));
// }

// ADReal
// LynxADLyakhovskyDamage::dYieldFunctiondMua()
// {
//   ADReal f = yieldFunction();

//   if (MooseUtils::absoluteFuzzyEqual(f, -1.0))
//     return 0.0;
//   else
//     return -_mua / (1.0 + f);
// }

// ADReal
// LynxADLyakhovskyDamage::d2YieldFunctiondDev2()
// {
//   ADReal f = yieldFunction();
//   ADReal df_ddev_stress = dYieldFunctiondDev();

//   if (MooseUtils::absoluteFuzzyEqual(_pressure, 0.0) || MooseUtils::absoluteFuzzyEqual(f, -1.0))
//     return 0.0;
//   else
//     return (1.0 / Utility::pow<2>(_mu0 * _pressure) - Utility::pow<2>(df_ddev_stress) / (1.0 + f)) /
//            (1.0 + f);
// }

// ADReal
// LynxADLyakhovskyDamage::d2YieldFunctiondDevdPres()
// {
//   ADReal f = yieldFunction();
//   ADReal df_ddev_stress = dYieldFunctiondDev();
//   ADReal df_dpressure = dYieldFunctiondPres();

//   if (MooseUtils::absoluteFuzzyEqual(_pressure, 0.0) || MooseUtils::absoluteFuzzyEqual(f, -1.0))
//     return 0.0;
//   else
//     return -df_ddev_stress * (2.0 / _pressure + df_dpressure / (1.0 + f));
// }

// ADReal
// LynxADLyakhovskyDamage::d2YieldFunctiondDevdMua()
// {
//   ADReal f = yieldFunction();
//   ADReal df_ddev_stress = dYieldFunctiondDev();
//   ADReal df_dmua = dYieldFunctiondMua();

//   if (MooseUtils::absoluteFuzzyEqual(f, -1.0))
//     return 0.0;
//   else
//     return -df_ddev_stress * df_dmua / (1.0 + f);
// }

// ADReal
// LynxADLyakhovskyDamage::d2YieldFunctiondDamdDev()
// {
//   ADReal f = yieldFunction();
//   ADReal df_ddamage_force = dYieldFunctiondDam();
//   ADReal df_ddev_stress = dYieldFunctiondDev();

//   if (MooseUtils::absoluteFuzzyEqual(f, -1.0))
//     return 0.0;
//   else
//     return -df_ddev_stress * df_ddamage_force / (1.0 + f);
// }

// ADReal
// LynxADLyakhovskyDamage::d2YieldFunctiondDamdPres()
// {
//   ADReal f = yieldFunction();
//   ADReal df_ddamage_force = dYieldFunctiondDam();
//   ADReal df_dpressure = dYieldFunctiondPres();

//   if (MooseUtils::absoluteFuzzyEqual(f, -1.0))
//     return 0.0;
//   else
//     return -df_dpressure * df_ddamage_force / (1.0 + f);
// }

// ADReal
// LynxADLyakhovskyDamage::d2YieldFunctiondDamdMua()
// {
//   ADReal f = yieldFunction();
//   ADReal df_ddamage_force = dYieldFunctiondDam();
//   ADReal df_dmua = dYieldFunctiondMua();

//   if (MooseUtils::absoluteFuzzyEqual(f, -1.0))
//     return 0.0;
//   else
//     return -(2.0 * _mua + df_dmua * df_ddamage_force / (1.0 + f)) / (1.0 + f);
// }

// ADReal
// LynxADLyakhovskyDamage::d2YieldFunctiondDam2()
// {
//   ADReal f = yieldFunction();
//   ADReal df_ddamage_force = dYieldFunctiondDam();

//   if (MooseUtils::absoluteFuzzyEqual(_damage_force, 0.0) || MooseUtils::absoluteFuzzyEqual(f, -1.0))
//     return 0.0;
//   else
//     return -df_ddamage_force * (1.0 / _damage_force + df_ddamage_force / (1.0 + f));
// }

// ADReal
// LynxADLyakhovskyDamage::strainRatio(const ADRankTwoTensor & elastic_strain)
// {
//   const ADReal strain_v = elastic_strain.trace();
//   const ADReal strain_norm = elastic_strain.L2norm();

//   if (strain_norm != 0.0)
//     return strain_v / strain_norm;
//   else
//     return -std::sqrt(3.0);
// }

// ADRankTwoTensor
// LynxADLyakhovskyDamage::reformPlasticStrainTensor(const ADReal & gamma_s)
// {
//   ADRankTwoTensor flow_dir =
//       (_dev_stress_tr != 0.0) ? _stress_tr.deviatoric() / _dev_stress_tr : ADRankTwoTensor();

//   ADRankTwoTensor delta_gamma = 1.5 * gamma_s * _dt * flow_dir;

//   return delta_gamma;
// }

// ADReal
// LynxADLyakhovskyDamage::computeDamagePoroMech(const ADReal & gamma_s, const ADReal & gamma_d, const ADRankTwoTensor & elastic_strain_incr)
// {
//   // Effective bulk modulus
//   ADReal Ke = _K - _damage_old[_qp] * _Gam0 * (2.0 / 3.0 * (_xi - _xi0) + _xi / 3.0  * (1.0 - Utility::pow<2>(_xi) / 3.0));

//   // Shear strain increment
//   ADReal dev_strain_rate = (_dt != 0.0) ? std::sqrt(2.0 / 3.0) * elastic_strain_incr.deviatoric().L2norm() / _dt : 0.0;

//   ADReal Gamv = (Ke != 0.0) ? _Gam0 / Ke * std::sqrt(1.5) * std::pow(1.0 - Utility::pow<2>(_xi) / 3.0, 1.5): 0.0;
//   ADReal Piv = (Ke != 0.0) ? _Gam0 / Ke * _e_norm * (1.0 + _xi / 3.0 * (_xi - 2.0 * _xi0)) : 0.0;

//   return _damage_old[_qp] * Gamv * (dev_strain_rate - gamma_s) + Piv * gamma_d;
// }