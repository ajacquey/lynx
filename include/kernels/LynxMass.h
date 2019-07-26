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

class LynxMass;

template <>
InputParameters validParams<LynxMass>();

class LynxMass : public DerivativeMaterialInterface<Kernel>
{
public:
  LynxMass(const InputParameters & parameters);
  virtual ~LynxMass() {}
  enum PenaltyType
  {
    LINEAR,
    LAPLACE
  };

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  const Real _penalty;
  PenaltyType _penalty_type;

  unsigned _ndisp;
  std::vector<unsigned> _disp_var;

  const MaterialProperty<RankTwoTensor> & _strain_increment;
};