/******************************************************************************/
/*                       LYNX, a MOOSE-based application                      */
/*                                                                            */
/*          Copyright (C) 2017 by Antoine B. Jacquey and Mauro Cacace         */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*    This program is free software: you can redistribute it and/or modify    */
/*    it under the terms of the GNU General Public License as published by    */
/*      the Free Software Foundation, either version 3 of the License, or     */
/*                     (at your option) any later version.                    */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program. If not, see <http://www.gnu.org/licenses/>     */
/******************************************************************************/

#ifndef LYNXSOLIDMOMENTUM_H
#define LYNXSOLIDMOMENTUM_H

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
  virtual void computeOffDiagJacobian(MooseVariableFEBase & jvar) override;
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

#endif // LYNXSOLIDMOMENTUM_H
