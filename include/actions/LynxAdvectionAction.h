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

#ifndef LYNXADVECTIONACTION_H
#define LYNXADVECTIONACTION_H

#include "Action.h"

class LynxAdvectionAction;

template <>
InputParameters validParams<LynxAdvectionAction>();

class LynxAdvectionAction : public Action
{
public:
  LynxAdvectionAction(InputParameters params);
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
  std::vector<VariableName> _displacements;
  std::vector<VariableName> _compositional_phases;
  std::vector<VariableName> _temperature;

  std::vector<std::string> _aux_variables;
  std::vector<std::string> _aux_kernels;
  std::vector<std::string> _max_var_pp;
  std::vector<std::string> _min_var_pp;
  std::vector<std::string> _avg_var_pp;
  std::vector<std::string> _max_entropy_pp;
  std::vector<std::string> _min_entropy_pp;
  std::vector<std::string> _avg_entropy_pp;
};

#endif // LYNXADVECTIONACTION_H
