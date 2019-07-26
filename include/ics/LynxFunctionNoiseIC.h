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

#ifndef LYNXFUNCTIONNOISEIC_H
#define LYNXFUNCTIONNOISEIC_H

#include "InitialCondition.h"

class LynxFunctionNoiseIC;
class Function;
class RandomIC;

template <>
InputParameters validParams<LynxFunctionNoiseIC>();

class LynxFunctionNoiseIC : public InitialCondition
{
public:
  LynxFunctionNoiseIC(const InputParameters & parameters);

protected:
  virtual Real value(const Point & p) override;
  virtual RealGradient gradient(const Point & p) override;

  const Function & _func;
  Real _rand_per;
};

#endif // LYNXFUNCTIONNOISEIC_H
