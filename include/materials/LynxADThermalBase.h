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

class LynxADThermalBase : public LynxADMaterialBase
{
public:
  static InputParameters validParams();
  LynxADThermalBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void computeQpSpecificHeat();
  virtual void computeQpThermalDiff();
  virtual void computeQpThermalExpansion();
  virtual void computeQpThermalSource();
  virtual void computeQpHeatCap() = 0;
  virtual void computeQpThermalCond() = 0;
  virtual void computeQpThermalExp() = 0;
  virtual Real computeMixtureProperty(const Real fluid_prop, const Real solid_prop);
  virtual ADReal computeADMixtureProperty(const ADReal fluid_prop, const ADReal solid_prop);

  // const bool _coupled_porosity;
  const VariableValue & _porosity;
  const std::vector<Real> _heat_source;
  const bool _coupled_dens;

  const ADMaterialProperty<Real> * _rho_f;
  const ADMaterialProperty<Real> * _rho_s;
  ADMaterialProperty<Real> & _thermal_diff;
  ADMaterialProperty<Real> & _rhoC_b;
  ADMaterialProperty<Real> & _rhoC_f;
  MaterialProperty<Real> & _thermal_exp;
  MaterialProperty<Real> & _radiogenic_heat;

  std::vector<Real> _c_f;
  std::vector<Real> _c_s;
  std::vector<Real> _lambda_f;
  std::vector<Real> _lambda_s;
  std::vector<Real> _beta_f;
  std::vector<Real> _beta_s;
};