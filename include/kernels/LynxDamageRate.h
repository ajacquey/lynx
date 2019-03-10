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
