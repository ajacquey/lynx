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
#include "LynxMass.h"

registerMooseObject("LynxApp", LynxMass);

template <>
InputParameters
validParams<LynxMass>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Divergence of solid velocity for incompressible Stoke flow.");
  params.addRequiredCoupledVar(
      "displacements", "The string of displacements variables suitable for the problem statement.");
  params.addParam<Real>("penalty", 0.0, "The value of the penalty.");
  MooseEnum penalty_type_options("linear=0 laplace=1", "linear");
  params.addParam<MooseEnum>(
      "penalty_type", penalty_type_options, "The type of penalty formulation.");

  return params;
}

LynxMass::LynxMass(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _penalty(getParam<Real>("penalty")),
    _penalty_type((PenaltyType)(int)parameters.get<MooseEnum>("penalty_type")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _strain_increment(getDefaultMaterialProperty<RankTwoTensor>("strain_increment"))
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

Real
LynxMass::computeQpResidual()
{
  Real one_on_penalty = (_penalty != 0) ? 1.0 / _penalty : 0.0;
  Real res = -_strain_increment[_qp].trace() * _test[_i][_qp];
  switch (_penalty_type)
  {
    case LINEAR:
      res += -one_on_penalty * _u[_qp] * _test[_i][_qp];
      break;
    case LAPLACE:
      res += -one_on_penalty * _grad_u[_qp] * _grad_test[_i][_qp];
      break;
  }

  return res;
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

Real
LynxMass::computeQpJacobian()
{
  Real one_on_penalty = (_penalty != 0) ? 1.0 / _penalty : 0.0;
  Real jac = 0.0;

  switch (_penalty_type)
  {
    case LINEAR:
      jac += -one_on_penalty * _phi[_j][_qp] * _test[_i][_qp];
      break;
    case LAPLACE:
      jac += -one_on_penalty * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
      break;
  }

  return jac;
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

Real
LynxMass::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real jac = 0.0;
  for (unsigned int coupled_component = 0; coupled_component < _ndisp; coupled_component++)
    if (jvar == _disp_var[coupled_component])
      jac += -_grad_phi[_j][_qp](coupled_component) * _test[_i][_qp] / _dt;

  return jac;
}
