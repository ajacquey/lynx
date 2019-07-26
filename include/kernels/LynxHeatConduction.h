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

#ifndef LYNXHEATCONDUCTION_H
#define LYNXHEATCONDUCTION_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

class LynxHeatConduction;

template <>
InputParameters validParams<LynxHeatConduction>();

class LynxHeatConduction : public DerivativeMaterialInterface<Kernel>
{
public:
  LynxHeatConduction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  // virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const MaterialProperty<Real> & _thermal_diff;
  const MaterialProperty<Real> & _rho;
  const MaterialProperty<Real> & _dinvrho_dtemp;
};

#endif // LYNXHEATCONDUCTION_H
