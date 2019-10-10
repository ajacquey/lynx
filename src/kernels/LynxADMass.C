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

#include "LynxADMass.h"

registerADMooseObject("LynxApp", LynxADMass);

defineADValidParams(
    LynxADMass,
    ADKernel,
    params.addClassDescription("Divergence of solid velocity for incompressible Stoke flow.");
    params.addParam<Real>("penalty", 0.0, "The value of the penalty.");
    MooseEnum penalty_type_options("linear=0 laplace=1", "linear");
    params.addParam<MooseEnum>("penalty_type",
                               penalty_type_options,
                               "The type of penalty formulation."););

template <ComputeStage compute_stage>
LynxADMass<compute_stage>::LynxADMass(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _penalty(getParam<Real>("penalty")),
    _penalty_type(getParam<MooseEnum>("penalty_type")),
    _strain_increment(getADMaterialProperty<RankTwoTensor>("strain_increment"))
{
}

template <ComputeStage compute_stage>
ADReal
LynxADMass<compute_stage>::computeQpResidual()
{
  Real one_on_penalty = (_penalty != 0.0) ? 1.0 / _penalty : 0.0;
  ADReal res = -_strain_increment[_qp].trace() * _test[_i][_qp];
  switch (_penalty_type)
  {
    case 0: // LINEAR
      res += -one_on_penalty * _u[_qp] * _test[_i][_qp];
      break;
    case 1: // LAPLACE
      res += -one_on_penalty * _grad_u[_qp] * _grad_test[_i][_qp];
      break;
  }

  return res;
}