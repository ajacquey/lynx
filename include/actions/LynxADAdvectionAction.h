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

class LynxADAdvectionAction;

template <>
InputParameters validParams<LynxADAdvectionAction>();

class LynxADAdvectionAction : public Action
{
public:
  LynxADAdvectionAction(InputParameters params);
  enum ElementLengthType
  {
    MIN,
    MAX,
    AVERAGE
  };
  enum ExecuteOnType
  {
    TIMESTEP_BEGIN,
    TIMESTEP_END
  };

  virtual void act() override;

protected:
  virtual void createAuxVariableActions();
  virtual void createAuxKernelActions();
  virtual void createPostProcessorActions();
  virtual void createKernelActions();

  ElementLengthType _element_length_type;
  ExecuteOnType _execute_on;

  Real _beta_stabilization;
  Real _cr_stabilization;
  Real _coeff_Hs;

  std::vector<VariableName> _velocities;
  std::vector<VariableName> _compositional_phases;
  std::vector<VariableName> _temperature;
  bool _has_inelastic_heat_var;
  std::vector<VariableName> _inelastic_heat;
  bool _has_pressure_var;
  std::vector<VariableName> _pressure;

  std::vector<std::string> _aux_variables;
  std::vector<std::string> _aux_kernels;
  std::vector<std::string> _max_var_pp;
  std::vector<std::string> _min_var_pp;
  std::vector<std::string> _avg_var_pp;
  std::vector<std::string> _max_entropy_pp;
  std::vector<std::string> _min_entropy_pp;
  std::vector<std::string> _avg_entropy_pp;
};