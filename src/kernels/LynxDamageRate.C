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

#include "LynxDamageRate.h"

registerMooseObject("LynxApp", LynxDamageRate);

template <>
InputParameters
validParams<LynxDamageRate>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(
      "Damage rate for damage rheology. Works fully coupled way or in a MultiApp");
  params.addCoupledVar("damage_rate", "The damage rate auxiliary variable");
  params.addCoupledVar("displacements", "The displacements variables.");
  return params;
}

LynxDamageRate::LynxDamageRate(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _u_old(valueOld()),
    _coupled_dam(isCoupled("damage_rate")),
    _coupled_disp(isCoupled("displacements")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _damage_rate(_coupled_dam ? coupledValue("damage_rate") : _zero),
    _damage_rate_mat(getDefaultMaterialProperty<Real>("damage_rate")),
    _ddamage_rate_dstrain(getDefaultMaterialProperty<RankTwoTensor>("ddamage_rate_dstrain"))
{
  if (_coupled_dam && _coupled_disp)
    mooseError("LynxDamageRate: you provided both damage_rate and displacements. Choose either "
               "damage_rate if running in a subApp or displacements if running a fully coupled "
               "simulation.");
  else if (!_coupled_dam && !_coupled_disp)
    mooseError("LynxDamageRate: you need to provide either damage_rate for running in a subApp or "
               "displacements for fully coupled!");

  for (unsigned i = 0; i < _ndisp; ++i)
    _disp_var[i] = _coupled_disp ? coupled("displacements", i) : -1;
}

/******************************************************************************/
/*                                  RESIDUALS                                 */
/******************************************************************************/

Real
LynxDamageRate::computeQpResidual()
{
  Real rate = (1.0 - _u_old[_qp]) / _dt;
  if (_coupled_dam)
    rate = std::min(_damage_rate[_qp], rate);
  else if (_coupled_disp)
    rate = std::min(_damage_rate_mat[_qp], rate);

  return -rate * _test[_i][_qp];
}

/******************************************************************************/
/*                                  JACOBIAN                                  */
/******************************************************************************/

Real
LynxDamageRate::computeQpJacobian()
{
  return 0.0;
}

/******************************************************************************/
/*                              OFF-DIAG JACOBIAN                             */
/******************************************************************************/

Real
LynxDamageRate::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real jac = 0.0;
  for (unsigned int coupled_component = 0; coupled_component < _ndisp; ++coupled_component)
    if (_coupled_disp && (jvar == _disp_var[coupled_component]))
    {
      Real rate = (1.0 - _u_old[_qp]) / _dt;
      if (_coupled_dam)
        rate = std::min(_damage_rate[_qp], rate);
      else if (_coupled_disp)
        rate = std::min(_damage_rate_mat[_qp], rate);

      if (rate != ((1.0 - _u_old[_qp]) / _dt))
        jac += -_ddamage_rate_dstrain[_qp].row(coupled_component) * _grad_phi[_j][_qp] *
               _test[_i][_qp];
    }

  return jac;
}
