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

#pragma once

#include "LynxADMaterialBase.h"

class LynxADHydroBase : public LynxADMaterialBase
{
public:
  static InputParameters validParams();
  LynxADHydroBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void computeQpCompressibilities();
  virtual void computeQpFluidMobility();
  virtual void computeQpPoroMech();
  virtual void computeQpFluidCompressibility() = 0;
  virtual void computeQpSolidCompressibility() = 0;
  virtual void computeQpPermeability() = 0;
  virtual void computeQpFluidViscosity() = 0;

  const VariableValue & _porosity;
  const bool _coupled_mech;
  const ADMaterialProperty<Real> * _K;
  const ADMaterialProperty<RankTwoTensor> * _strain_increment;
  const bool _has_viscous;
  const ADMaterialProperty<RankTwoTensor> * _viscous_strain_incr;
  const bool _has_plastic;
  const ADMaterialProperty<RankTwoTensor> * _plastic_strain_incr;
  ADMaterialProperty<Real> & _biot;
  ADMaterialProperty<Real> & _C_d;
  ADMaterialProperty<Real> & _C_biot;
  ADMaterialProperty<Real> & _fluid_mobility;
  ADMaterialProperty<Real> & _poro_mech;

  std::vector<Real> _C_f;
  std::vector<Real> _C_s;
  std::vector<Real> _k;
  std::vector<Real> _eta_f;
};