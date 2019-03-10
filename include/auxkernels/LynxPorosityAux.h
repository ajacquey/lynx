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

#ifndef LYNXPOROSITYAUX_H
#define LYNXPOROSITYAUX_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "DerivativeMaterialInterface.h"

class LynxPorosityAux;

template <>
InputParameters validParams<LynxPorosityAux>();

class LynxPorosityAux : public DerivativeMaterialInterface<AuxKernel>
{
public:
  LynxPorosityAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  virtual Real computeEvDot();
  virtual Real computeEvInDot();

  bool _coupled_pf;
  const VariableValue & _pf_dot;
  const MaterialProperty<Real> & _biot;
  const MaterialProperty<Real> & _C_d;
  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _inelastic_strain;
  const MaterialProperty<RankTwoTensor> & _inelastic_strain_old;
};

#endif // LYNXPOROSITYAUX_H
