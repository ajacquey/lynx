//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

namespace libMesh
{
template <typename>
class VectorValue;
typedef VectorValue<Real> RealGradient;
}
using libMesh::RealGradient;

namespace ElasticityTensorTools
{
/**
 * Get the shear modulus for an isotropic elasticity tensor
 * param elasticity_tensor the tensor (must be isotropic, but not checked for efficiency)
 */
template <typename T>
T
getIsotropicShearModulus(const RankFourTensorTempl<T> & elasticity_tensor)
{
  return elasticity_tensor(0, 1, 0, 1);
}

/**
 * Get the bulk modulus for an isotropic elasticity tensor
 * param elasticity_tensor the tensor (must be isotropic, but not checked for efficiency)
 */
template <typename T>
T
getIsotropicBulkModulus(const RankFourTensorTempl<T> & elasticity_tensor)
{
  const T shear_modulus = getIsotropicShearModulus(elasticity_tensor);
  // dilatational modulus is defined as lambda plus two mu
  const T dilatational_modulus = elasticity_tensor(0, 0, 0, 0);
  const T lambda = dilatational_modulus - 2.0 * shear_modulus;
  const T bulk_modulus = lambda + 2.0 * shear_modulus / 3.0;
  return bulk_modulus;
}

/**
 * Get the Young's modulus for an isotropic elasticity tensor
 * param elasticity_tensor the tensor (must be isotropic, but not checked for efficiency)
 */
template <typename T>
T
getIsotropicYoungsModulus(const RankFourTensorTempl<T> & elasticity_tensor)
{
  const T shear_modulus = getIsotropicShearModulus(elasticity_tensor);
  // dilatational modulus is defined as lambda plus two mu
  const T dilatational_modulus = elasticity_tensor(0, 0, 0, 0);
  const T lambda = dilatational_modulus - 2.0 * shear_modulus;
  const T youngs_modulus =
      shear_modulus * (3.0 * lambda + 2.0 * shear_modulus) / (lambda + shear_modulus);
  return youngs_modulus;
}

/**
 * Get the Poisson's modulus for an isotropic elasticity tensor
 * param elasticity_tensor the tensor (must be isotropic, but not checked for efficiency)
 */
template <typename T>
T
getIsotropicPoissonsRatio(const RankFourTensorTempl<T> & elasticity_tensor)
{
  const T poissons_ratio = elasticity_tensor(1, 1, 0, 0) /
                           (elasticity_tensor(1, 1, 1, 1) + elasticity_tensor(1, 1, 0, 0));
  return poissons_ratio;
}

/**
 * Construct the elasticity tensor based on the bulk and shear modulus
 */
template <typename T>
RankFourTensorTempl<T>
elasticityTensorKandG(const T & K, const T & G)
{
  std::vector<T> iso_const(2);
  iso_const[0] = K - 2.0 / 3.0 * G;
  iso_const[1] = G;

  return RankFourTensorTempl<T>(iso_const, RankFourTensorTempl<T>::symmetric_isotropic);
}

}