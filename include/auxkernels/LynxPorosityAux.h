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