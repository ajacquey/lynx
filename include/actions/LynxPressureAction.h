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

#include "Action.h"

class LynxPressureAction;

template <>
InputParameters validParams<LynxPressureAction>();

class LynxPressureAction : public Action
{
public:
  LynxPressureAction(const InputParameters & params);

  virtual void act() override;

protected:
  std::vector<std::vector<AuxVariableName>> _save_in_vars;
  std::vector<bool> _has_save_in_vars;
};