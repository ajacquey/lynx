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

#ifndef LYNXHOLDSTRESSBC_H
#define LYNXHOLDSTRESSBC_H

#include "IntegratedBC.h"

class LynxHoldStressBC;
class Function;

template <>
InputParameters validParams<LynxHoldStressBC>();

class LynxHoldStressBC : public IntegratedBC
{
public:
  LynxHoldStressBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  bool _coupled_pf;
  const VariableValue & _pf;
  const int _component;
  const Real _bulk_modulus;
  const Real _shear_modulus;
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  const MaterialProperty<Real> & _biot_coeff;
};

#endif // LYNXHOLDSTRESSBC_H
