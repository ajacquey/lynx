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

#ifndef LYNXDAMAGERATE_H
#define LYNXDAMAGERATE_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

class LynxDamageRate;

template <>
InputParameters validParams<LynxDamageRate>();

class LynxDamageRate : public DerivativeMaterialInterface<Kernel>
{
public:
  LynxDamageRate(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _u_old;
  bool _coupled_dam;
  bool _coupled_disp;
  unsigned int _ndisp;
  std::vector<unsigned> _disp_var;
  const VariableValue & _damage_rate;
  const MaterialProperty<Real> & _damage_rate_mat;
  const MaterialProperty<RankTwoTensor> & _ddamage_rate_dstrain;
};

#endif // LYNXDAMAGERATE_H
