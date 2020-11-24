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

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class LynxSolidMomentum;

template <>
InputParameters validParams<LynxSolidMomentum>();

class LynxSolidMomentum : public DerivativeMaterialInterface<Kernel>
{
public:
  LynxSolidMomentum(const InputParameters & parameters);

  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(unsigned int jvar) override;
  using Kernel::computeOffDiagJacobian;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  virtual Real elasticJacobian(const RankFourTensor & jacobian_r4t,
                               unsigned int i,
                               unsigned int k,
                               const RealGradient & grad_test,
                               const RealGradient & grad_phi);
  virtual void computeAverageGradientTest();
  virtual void computeAverageGradientPhi();

  unsigned int _ndisp;
  bool _coupled_temp;
  bool _coupled_pf;
  const VariableValue & _pf;
  bool _coupled_plith;
  bool _coupled_pdyn;
  bool _coupled_dam;
  const unsigned int _component;
  bool _vol_locking_correction;
  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<Real> & _biot;
  const MaterialProperty<RankFourTensor> & _tangent_modulus;
  const MaterialProperty<RankTwoTensor> & _dthermal_strain_dtemp;
  const MaterialProperty<RealVectorValue> & _gravity;
  const MaterialProperty<Real> & _rho_b;
  const MaterialProperty<Real> & _drho_dtemp;
  const MaterialProperty<Real> & _drho_dev;
  const MaterialProperty<RankTwoTensor> & _dstress_ddamage;
  // Gradient of test function averaged over the element. Used in volumetric locking correction
  // calculation.
  std::vector<std::vector<Real>> _avg_grad_test;
  // Gradient of phi function averaged over the element. Used in volumetric locking correction
  // calculation.
  std::vector<std::vector<Real>> _avg_grad_phi;
  std::vector<unsigned int> _disp_var;
  unsigned int _temp_var;
  unsigned int _pf_var;
  unsigned int _plith_var;
  unsigned int _pdyn_var;
  unsigned int _damage_var;
};