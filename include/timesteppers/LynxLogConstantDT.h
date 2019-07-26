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

#ifndef LYNXLOGCONSTANTDT_H
#define LYNXLOGCONSTANTDT_H

#include "TimeStepper.h"

class LynxLogConstantDT;

template <>
InputParameters validParams<LynxLogConstantDT>();

class LynxLogConstantDT : public TimeStepper
{
public:
  LynxLogConstantDT(const InputParameters & parameters);

protected:
  virtual Real computeInitialDT() override;
  virtual Real computeDT() override;

private:
  const Real _log_dt;
  const Real _first_dt;
  const Real _max_dt;
  const Real _dt_factor;
  const Real _growth_factor;
};

#endif // LYNXLOGCONSTANTDT_H
