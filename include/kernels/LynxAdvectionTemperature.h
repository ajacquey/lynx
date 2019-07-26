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

#ifndef LYNXADVECTIONTEMPERATURE_H
#define LYNXADVECTIONTEMPERATURE_H

#include "LynxAdvectionBase.h"

class LynxAdvectionTemperature;

template <>
InputParameters validParams<LynxAdvectionTemperature>();

class LynxAdvectionTemperature : public LynxAdvectionBase
{
public:
  LynxAdvectionTemperature(const InputParameters & parameters);

protected:
  virtual Real computeArtificialViscosity() override;
  virtual void computeEntropyResidual();

  // diffusion residual
  const MaterialProperty<Real> & _thermal_diff;
  // source term residual
  Real _coeff_Hs;
  const MaterialProperty<Real> & _rhoC;
  const MaterialProperty<Real> & _radiogenic_heat;
  const MaterialProperty<Real> & _inelastic_heat;
};

#endif // LYNXADVECTIONTEMPERATURE_H
