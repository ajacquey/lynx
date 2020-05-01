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

#include "ADKernel.h"

class LynxADSolidMomentum : public ADKernel
{
public:
  static InputParameters validParams();
  LynxADSolidMomentum(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
  void precalculateResidual() override;

  const ADVariableValue & _pf;
  const unsigned int _component;
  const bool _vol_locking_correction;
  const ADMaterialProperty<RankTwoTensor> & _stress;
  const bool _coupled_pf;
  const ADMaterialProperty<Real> * _biot;
  const bool _coupled_grav;
  const ADMaterialProperty<RealVectorValue> * _gravity;
  const ADMaterialProperty<Real> * _rho_b;
  std::vector<ADReal> _avg_grad_test;
};