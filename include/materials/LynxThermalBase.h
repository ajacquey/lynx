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

#include "LynxMaterialBase.h"
#include "RankTwoTensor.h"

class LynxThermalBase;

template <>
InputParameters validParams<LynxThermalBase>();

class LynxThermalBase : public LynxMaterialBase
{
public:
  LynxThermalBase(const InputParameters & parameters);
  virtual ~LynxThermalBase() {}

protected:
  virtual void computeQpProperties() override;
  virtual void computeQpSpecificHeat();
  virtual void computeQpThermalDiff();
  virtual void computeQpThermalExpansion();
  virtual void computeQpThermalStrain();
  virtual void computeQpThermalSource();
  virtual void computeQpHeatCap() = 0;
  virtual void computeQpThermalCond() = 0;
  virtual void computeQpThermalExp() = 0;
  virtual Real computeMixtureProperty(const Real fluid_prop, const Real solid_prop);

  const VariableValue & _temp_dot;
  const VariableValue & _porosity;
  const std::vector<Real> _heat_source;

  const MaterialProperty<Real> & _rho_f;
  const MaterialProperty<Real> & _rho_s;
  MaterialProperty<Real> & _thermal_diff;
  MaterialProperty<Real> & _rhoC_b;
  MaterialProperty<Real> & _rhoC_f;
  MaterialProperty<Real> & _thermal_exp;
  MaterialProperty<RankTwoTensor> & _thermal_strain_incr;
  MaterialProperty<RankTwoTensor> & _dthermal_strain_dtemp;
  MaterialProperty<Real> & _radiogenic_heat;

  std::vector<Real> _c_f;
  std::vector<Real> _c_s;
  std::vector<Real> _lambda_f;
  std::vector<Real> _lambda_s;
  std::vector<Real> _beta_f;
  std::vector<Real> _beta_s;
};