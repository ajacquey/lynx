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

#ifndef LYNXHEATFLUXAUX_H
#define LYNXHEATFLUXAUX_H

#include "AuxKernel.h"
#include "DerivativeMaterialInterface.h"

class LynxHeatFluxAux;

template <>
InputParameters validParams<LynxHeatFluxAux>();

class LynxHeatFluxAux : public DerivativeMaterialInterface<AuxKernel>
{
public:
  LynxHeatFluxAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const unsigned int _component;
  const VariableGradient & _grad_T;
  const MaterialProperty<Real> & _rhoC_b;
  const MaterialProperty<Real> & _thermal_diff;
  Real _cond_bulk;
};

#endif // LynxHeatFluxAux_H
