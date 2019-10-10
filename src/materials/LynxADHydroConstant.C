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

#include "LynxADHydroConstant.h"

registerADMooseObject("LynxApp", LynxADHydroConstant);

defineADValidParams(
    LynxADHydroConstant,
    LynxADHydroBase,
    params.addClassDescription("Constant hydraulic properties.");
    params.addRequiredParam<std::vector<Real>>("permeability", "The permeability of the matrix.");
    params.addRequiredParam<std::vector<Real>>("fluid_viscosity", "The viscosity of the fluid.");
    params.addParam<std::vector<Real>>("fluid_modulus", "The bulk modulus of the fluid phase.");
    params.addParam<std::vector<Real>>("solid_modulus", "The bulk modulus of the solid phase."););

template <ComputeStage compute_stage>
LynxADHydroConstant<compute_stage>::LynxADHydroConstant(const InputParameters & parameters)
  : LynxADHydroBase<compute_stage>(parameters),
    _perm(this->getLynxParam("permeability")),
    _fluid_viscosity(this->getLynxParam("fluid_viscosity"))
{
  if (isParamValid("fluid_modulus"))
  {
    _fluid_compr = std::vector<Real>(_n_composition, 0.0);
    std::vector<Real> fluid_modulus = this->getLynxParam("fluid_modulus");
    for (unsigned int i = 0; i < _n_composition; ++i)
      if (fluid_modulus[i] != 0.0)
        _fluid_compr[i] = 1.0 / fluid_modulus[i];
  }
  else
    _fluid_compr = std::vector<Real>(_n_composition, 0.0);

  if (isParamValid("solid_modulus"))
  {
    _solid_compr = std::vector<Real>(_n_composition, 0.0);
    std::vector<Real> solid_modulus = this->getLynxParam("solid_modulus");
    for (unsigned int i = 0; i < _n_composition; ++i)
      if (solid_modulus[i] != 0.0)
        _solid_compr[i] = 1.0 / solid_modulus[i];
  }
  else
    _solid_compr = std::vector<Real>(_n_composition, 0.0);
}

template <ComputeStage compute_stage>
void
LynxADHydroConstant<compute_stage>::computeQpFluidCompressibility()
{
  _C_f[_qp] = this->averageProperty(_fluid_compr);
}

template <ComputeStage compute_stage>
void
LynxADHydroConstant<compute_stage>::computeQpSolidCompressibility()
{
  _C_s[_qp] = this->averageProperty(_solid_compr);
}


template <ComputeStage compute_stage>
void
LynxADHydroConstant<compute_stage>::computeQpPermeability()
{
  _k[_qp] = this->averageProperty(_perm);
}


template <ComputeStage compute_stage>
void
LynxADHydroConstant<compute_stage>::computeQpFluidViscosity()
{
  _eta_f[_qp] = this->averageProperty(_fluid_viscosity);
}

adBaseClass(LynxADHydroConstant);