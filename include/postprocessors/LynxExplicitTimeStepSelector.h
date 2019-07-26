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
#ifndef LYNXEXPLICITTIMESTEPSELECTOR_H
#define LYNXEXPLICITTIMESTEPSELECTOR_H

#include "ElementPostprocessor.h"
#include "DerivativeMaterialInterface.h"

class LynxExplicitTimeStepSelector;

template <>
InputParameters validParams<LynxExplicitTimeStepSelector>();

class LynxExplicitTimeStepSelector : public DerivativeMaterialInterface<ElementPostprocessor>
{
public:
  LynxExplicitTimeStepSelector(const InputParameters & parameters);
  virtual ~LynxExplicitTimeStepSelector();

  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
  virtual void threadJoin(const UserObject & uo);

protected:
  Real _value;
  const VariableValue & _vel_norm;
  Real _beta;
  Real _epsilon;
  bool _has_premult;
  Real _initial_value;
  Real _maximum_value;
};

#endif // LYNXEXPLICITTIMESTEPSELECTOR_H
