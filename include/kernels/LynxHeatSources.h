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

#ifndef LYNXHEATSOURCES_H
#define LYNXHEATSOURCES_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"
#include "RankTwoTensor.h"

class LynxHeatSources;

template <>
InputParameters validParams<LynxHeatSources>();

class LynxHeatSources : public DerivativeMaterialInterface<Kernel>
{
public:
  LynxHeatSources(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  Real _coeff_Hs;

  const MaterialProperty<Real> & _rhoC_b;
  const MaterialProperty<Real> & _radiogenic_heat;
  const MaterialProperty<Real> & _rho_b;
  const MaterialProperty<Real> & _dinvrho_dtemp;
  const MaterialProperty<Real> & _inelastic_heat_mat;
  const MaterialProperty<Real> & _adiabatic_heat;
  const MaterialProperty<Real> & _damage_heat;
  // const MaterialProperty<Real> & _dinelastic_heat_dtemp;
  const MaterialProperty<RankTwoTensor> & _dinelastic_heat_dstrain;
  const bool _coupled_disp;
  unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;
  const bool _coupled_inelastic_heat;
  const VariableValue & _inelastic_heat;
};

#endif // LYNXHEATSOURCES_H
