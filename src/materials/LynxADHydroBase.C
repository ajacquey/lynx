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

#include "LynxADHydroBase.h"

defineADValidParams(
    LynxADHydroBase,
    LynxADMaterialBase,
    params.addClassDescription("Base class for calculating the thermal properties.");
    params.addRequiredCoupledVar("porosity", "The porosity auxiliary variable."););

template <ComputeStage compute_stage>
LynxADHydroBase<compute_stage>::LynxADHydroBase(const InputParameters & parameters)
  : LynxADMaterialBase<compute_stage>(parameters),
    _porosity(adCoupledValue("porosity")),
    _coupled_mech(hasMaterialProperty<Real>("bulk_modulus")),
    _K(_coupled_mech ? &getADMaterialProperty<Real>("bulk_modulus") : nullptr),
    _strain_increment(_coupled_mech ? &getADMaterialProperty<RankTwoTensor>("strain_increment") : nullptr),
    // _viscous_strain_incr(_coupled_mech ? &getADMaterialProperty<RankTwoTensor>("viscous_strain_increment") : nullptr),
    // _plastic_strain_incr(_coupled_mech ? &getADMaterialProperty<RankTwoTensor>("plastic_strain_increment") : nullptr),
    _biot(declareADProperty<Real>("biot_coefficient")),
    _C_d(declareADProperty<Real>("bulk_compressibility")),
    _C_biot(declareADProperty<Real>("biot_compressibility")),
    _fluid_mobility(declareADProperty<Real>("fluid_mobility")),
    _poro_mech(declareADProperty<Real>("poro_mechanical")),
    _C_f(_fe_problem.getMaxQps()),
    _C_s(_fe_problem.getMaxQps()),
    _k(_fe_problem.getMaxQps()),
    _eta_f(_fe_problem.getMaxQps())
{
}

template <ComputeStage compute_stage>
void
LynxADHydroBase<compute_stage>::computeQpProperties()
{
  computeQpCompressibilities();
  computeQpFluidMobility();
  computeQpPoroMech();
}

template <ComputeStage compute_stage>
void
LynxADHydroBase<compute_stage>::computeQpCompressibilities()
{
  computeQpFluidCompressibility();
  computeQpSolidCompressibility();

  // Drained compressibility
  _C_d[_qp] = (_coupled_mech && ((*_K)[_qp] != 0.0)) ? 1.0 / (*_K)[_qp] : 0.0;

  // Biot coefficient
  _biot[_qp] = 1.0;
  if (_coupled_mech && _C_d[_qp] != 0.0)
    _biot[_qp] -= _C_s[_qp] / _C_d[_qp];

  // Pore compressibility
  ADReal C_phi = (_biot[_qp] - _porosity[_qp]) * _C_d[_qp];

  // Biot compressibility
  _C_biot[_qp] = _porosity[_qp] * _C_f[_qp] + (1.0 - _biot[_qp]) * C_phi;
}

template <ComputeStage compute_stage>
void
LynxADHydroBase<compute_stage>::computeQpFluidMobility()
{
  computeQpPermeability();
  computeQpFluidViscosity();

  // Fluid mobility
  _fluid_mobility[_qp] = _k[_qp] / _eta_f[_qp];
  if (_C_biot[_qp] != 0.0)
    _fluid_mobility[_qp] /= _C_biot[_qp];
}

template <ComputeStage compute_stage>
void
LynxADHydroBase<compute_stage>::computeQpPoroMech()
{
  if (_coupled_mech)
  {
    ADRankTwoTensor e_tot = (*_strain_increment)[_qp] / _dt;
    // ADRankTwoTensor e_in = ((*_viscous_strain_incr)[_qp] + (*_plastic_strain_incr)[_qp]) / _dt;
    ADRankTwoTensor e_in = ADRankTwoTensor();
    _poro_mech[_qp] = _biot[_qp] * e_tot.trace() + (1.0 - _biot[_qp]) * e_in.trace();

    if (_C_biot[_qp] != 0.0)
      _poro_mech[_qp] /= _C_biot[_qp];
  }
  else
    _poro_mech[_qp] = 0.0;
}

adBaseClass(LynxADHydroBase);