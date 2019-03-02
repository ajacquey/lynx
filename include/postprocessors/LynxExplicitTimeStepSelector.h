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
