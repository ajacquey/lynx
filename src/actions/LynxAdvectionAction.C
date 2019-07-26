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


#include "LynxAdvectionAction.h"
#include "Factory.h"
#include "FEProblem.h"
#include "MooseMesh.h"

#include "libmesh/fe_type.h"
#include "libmesh/string_to_enum.h"

registerMooseAction("LynxApp", LynxAdvectionAction, "add_aux_variable");
registerMooseAction("LynxApp", LynxAdvectionAction, "add_aux_kernel");
registerMooseAction("LynxApp", LynxAdvectionAction, "add_postprocessor");
registerMooseAction("LynxApp", LynxAdvectionAction, "add_kernel");

template <>
InputParameters
validParams<LynxAdvectionAction>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<std::vector<VariableName>>(
      "velocities", "The name of the advective velocities in the simulation.");
  params.addParam<std::vector<VariableName>>(
      "displacements", "The name of the advective displacements in the simulation.");
  params.addParam<std::vector<VariableName>>(
      "compositional_phases", "The name of the compositional phases to be advected.");
  params.addParam<std::vector<VariableName>>(
      "temperature", "The name of the temperature variable to be advected.");
  MooseEnum element_length_type_options("min=0 max=1 average=2", "min");
  params.addParam<MooseEnum>(
      "element_length_type", element_length_type_options, "The diameter of a single cell.");
  MooseEnum execute_on_options("timestep_begin=0 timestep_end=1", "timestep_begin");
  params.addParam<MooseEnum>("execute_on",
                             execute_on_options,
                             "The time when the postprocessors are asked to be executed.");
  params.addParam<Real>("beta_stabilization", 0.026, "The beta local stabilization parameter.");
  params.addParam<Real>("cr_stabilization", 0.5, "The cr local stabilization parameter.");
  params.addParam<Real>(
      "coeff_shear_heating", 0.0, "The coefficient in front of the shear heating generation.");
  return params;
}

LynxAdvectionAction::LynxAdvectionAction(InputParameters params)
  : Action(params),
    _element_length_type((ElementLengthType)(int)params.get<MooseEnum>("element_length_type")),
    _execute_on((ExecuteOnType)(int)params.get<MooseEnum>("execute_on")),
    _beta_stabilization(getParam<Real>("beta_stabilization")),
    _cr_stabilization(getParam<Real>("cr_stabilization")),
    _coeff_Hs(getParam<Real>("coeff_shear_heating")),
    _velocities(getParam<std::vector<VariableName>>("velocities")),
    _displacements(getParam<std::vector<VariableName>>("displacements"))
{
  if (isParamValid("compositional_phases"))
    _compositional_phases = getParam<std::vector<VariableName>>("compositional_phases");
  if (isParamValid("temperature"))
    _temperature = getParam<std::vector<VariableName>>("temperature");
  if (!isParamValid("temperature") && !isParamValid("compositional_phases"))
    mooseError("LynxAdvectionAction: no advected quantities are set.");
}

void
LynxAdvectionAction::act()
{
  if (_current_task == "add_aux_variable")
    createAuxVariableActions();
  else if (_current_task == "add_postprocessor")
    createPostProcessorActions();
  else if (_current_task == "add_aux_kernel")
    createAuxKernelActions();
  else if (_current_task == "add_kernel")
    createKernelActions();
}

/******************************************************************************/
/*                          Add Aux Variable                                   */
/******************************************************************************/

void
LynxAdvectionAction::createAuxVariableActions()
{
  const bool second_order = _problem->mesh().hasSecondOrderElements();
  for (unsigned i = 0; i < _compositional_phases.size(); ++i)
  {
    _aux_variables.push_back("entropy_" + _compositional_phases[i]);
    _problem->addAuxVariable(
        _aux_variables.back(),
        FEType(Utility::string_to_enum<Order>(second_order ? "SECOND" : "FIRST"),
               Utility::string_to_enum<FEFamily>("LAGRANGE")));
  }
  for (unsigned i = 0; i < _temperature.size(); ++i)
  {
    _aux_variables.push_back("entropy_" + _temperature[i]);
    _problem->addAuxVariable(
        _aux_variables.back(),
        FEType(Utility::string_to_enum<Order>(second_order ? "SECOND" : "FIRST"),
               Utility::string_to_enum<FEFamily>("LAGRANGE")));
  }
}

/******************************************************************************/
/*                          Add PostProcessor                                 */
/******************************************************************************/
void
LynxAdvectionAction::createPostProcessorActions()
{
  std::vector<VariableName> tmp; // this sucks, but MOOSE requires a std::vector each time
  tmp.resize(1);
  InputParameters params = emptyInputParameters();
  // velocity Postprocessors
  {
    params = _factory.getValidParams("LynxExtremeVectorValue");
    tmp[0] = _velocities[0];
    params.set<std::vector<VariableName>>("variable") = tmp;
    if (_velocities.size() > 1)
    {
      tmp[0] = _velocities[1];
      params.set<std::vector<VariableName>>("add_var_1") = tmp;
    }
    if (_velocities.size() > 2)
    {
      tmp[0] = _velocities[2];
      params.set<std::vector<VariableName>>("add_var_2") = tmp;
    }
    params.set<MooseEnum>("value_type") = "MAX";
    switch (_execute_on)
    {
      case TIMESTEP_BEGIN:
        params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
        break;
      case TIMESTEP_END:
        params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
        break;
    }
    _problem->addPostprocessor("LynxExtremeVectorValue", "max_vel_pp", params);
  }
  // compositional phases Postprocessors
  for (unsigned i = 0; i < _compositional_phases.size(); ++i)
  {
    _max_var_pp.push_back("max_" + _compositional_phases[i] + "_pp");
    _min_var_pp.push_back("min_" + _compositional_phases[i] + "_pp");
    _avg_var_pp.push_back("avg_" + _compositional_phases[i] + "_pp");
    _max_entropy_pp.push_back("max_entropy_" + _compositional_phases[i] + "_pp");
    _min_entropy_pp.push_back("min_entropy_" + _compositional_phases[i] + "_pp");
    _avg_entropy_pp.push_back("avg_entropy_" + _compositional_phases[i] + "_pp");
    tmp[0] = _compositional_phases[i];
    // max phase
    {
      params = _factory.getValidParams("LynxExtrapolateValue");
      params.set<std::vector<VariableName>>("variable") = tmp;
      switch (_execute_on)
      {
        case TIMESTEP_BEGIN:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
          break;
        case TIMESTEP_END:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          break;
      }
      params.set<MooseEnum>("value_type") = "MAX";
      _problem->addPostprocessor("LynxExtrapolateValue", _max_var_pp.back(), params);
    }
    // min phase
    {
      params = _factory.getValidParams("LynxExtrapolateValue");
      params.set<std::vector<VariableName>>("variable") = tmp;
      switch (_execute_on)
      {
        case TIMESTEP_BEGIN:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
          break;
        case TIMESTEP_END:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          break;
      }
      params.set<MooseEnum>("value_type") = "MIN";
      _problem->addPostprocessor("LynxExtrapolateValue", _min_var_pp.back(), params);
    }
    // average phase
    {
      params = _factory.getValidParams("LynxElementAverageValue");
      params.set<std::vector<VariableName>>("variable") = tmp;
      switch (_execute_on)
      {
        case TIMESTEP_BEGIN:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
          break;
        case TIMESTEP_END:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          break;
      }
      _problem->addPostprocessor("LynxElementAverageValue", _avg_var_pp.back(), params);
    }
    tmp[0] = _aux_variables[i];
    // max entropy phase
    {
      params = _factory.getValidParams("LynxExtrapolateValue");
      params.set<std::vector<VariableName>>("variable") = tmp;
      switch (_execute_on)
      {
        case TIMESTEP_BEGIN:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
          break;
        case TIMESTEP_END:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          break;
      }
      params.set<MooseEnum>("value_type") = "MAX";
      _problem->addPostprocessor("LynxExtrapolateValue", _max_entropy_pp.back(), params);
    }
    // min entropy phase
    {
      params = _factory.getValidParams("LynxExtrapolateValue");
      params.set<std::vector<VariableName>>("variable") = tmp;
      switch (_execute_on)
      {
        case TIMESTEP_BEGIN:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
          break;
        case TIMESTEP_END:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          break;
      }
      params.set<MooseEnum>("value_type") = "MIN";
      _problem->addPostprocessor("LynxExtrapolateValue", _min_entropy_pp.back(), params);
    }
    // average entropy phase
    {
      params = _factory.getValidParams("LynxElementAverageValue");
      params.set<std::vector<VariableName>>("variable") = tmp;
      switch (_execute_on)
      {
        case TIMESTEP_BEGIN:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
          break;
        case TIMESTEP_END:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          break;
      }
      _problem->addPostprocessor("LynxElementAverageValue", _avg_entropy_pp.back(), params);
    }
  }
  // temperature Postprocessors
  for (unsigned i = 0; i < _temperature.size(); ++i)
  {
    _max_var_pp.push_back("max_" + _temperature[i] + "_pp");
    _min_var_pp.push_back("min_" + _temperature[i] + "_pp");
    _avg_var_pp.push_back("avg_" + _temperature[i] + "_pp");
    _max_entropy_pp.push_back("max_entropy_" + _temperature[i] + "_pp");
    _min_entropy_pp.push_back("min_entropy_" + _temperature[i] + "_pp");
    _avg_entropy_pp.push_back("avg_entropy_" + _temperature[i] + "_pp");
    tmp[0] = _temperature[i];
    // max temperature
    {
      params = _factory.getValidParams("LynxExtrapolateValue");
      params.set<std::vector<VariableName>>("variable") = tmp;
      switch (_execute_on)
      {
        case TIMESTEP_BEGIN:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
          break;
        case TIMESTEP_END:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          break;
      }
      params.set<MooseEnum>("value_type") = "MAX";
      _problem->addPostprocessor("LynxExtrapolateValue", _max_var_pp.back(), params);
    }
    // min temperature
    {
      params = _factory.getValidParams("LynxExtrapolateValue");
      params.set<std::vector<VariableName>>("variable") = tmp;
      switch (_execute_on)
      {
        case TIMESTEP_BEGIN:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
          break;
        case TIMESTEP_END:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          break;
      }
      params.set<MooseEnum>("value_type") = "MIN";
      _problem->addPostprocessor("LynxExtrapolateValue", _min_var_pp.back(), params);
    }
    // average temperature
    {
      params = _factory.getValidParams("LynxElementAverageValue");
      params.set<std::vector<VariableName>>("variable") = tmp;
      switch (_execute_on)
      {
        case TIMESTEP_BEGIN:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
          break;
        case TIMESTEP_END:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          break;
      }
      _problem->addPostprocessor("ElementAverageValue", _avg_var_pp.back(), params);
    }
    tmp[0] = _aux_variables.back();
    // max entropy temperature
    {
      params = _factory.getValidParams("LynxExtrapolateValue");
      params.set<std::vector<VariableName>>("variable") = tmp;
      switch (_execute_on)
      {
        case TIMESTEP_BEGIN:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
          break;
        case TIMESTEP_END:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          break;
      }
      params.set<MooseEnum>("value_type") = "MAX";
      _problem->addPostprocessor("LynxExtrapolateValue", _max_entropy_pp.back(), params);
    }
    // min entropy temperature
    {
      params = _factory.getValidParams("LynxExtrapolateValue");
      params.set<std::vector<VariableName>>("variable") = tmp;
      switch (_execute_on)
      {
        case TIMESTEP_BEGIN:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
          break;
        case TIMESTEP_END:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          break;
      }
      params.set<MooseEnum>("value_type") = "MIN";
      _problem->addPostprocessor("LynxExtrapolateValue", _min_entropy_pp.back(), params);
    }
    // average entropy temperature
    {
      params = _factory.getValidParams("LynxElementAverageValue");
      params.set<std::vector<VariableName>>("variable") = tmp;
      switch (_execute_on)
      {
        case TIMESTEP_BEGIN:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
          break;
        case TIMESTEP_END:
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          break;
      }
      _problem->addPostprocessor("LynxElementAverageValue", _avg_entropy_pp.back(), params);
    }
  }
}

/******************************************************************************/
/*                          Add Aux Kernels                                   */
/******************************************************************************/
void
LynxAdvectionAction::createAuxKernelActions()
{
  std::vector<VariableName> tmp;
  tmp.resize(1);
  std::string type = "LynxEntropyAux";
  InputParameters params = emptyInputParameters();
  for (unsigned i = 0; i < _compositional_phases.size(); ++i)
  {
    tmp[0] = _compositional_phases[i];
    _aux_kernels.push_back("entropy_" + _compositional_phases[i] + "_aux");
    params = _factory.getValidParams(type);
    params.set<AuxVariableName>("variable") = _aux_variables[i];
    params.set<std::vector<VariableName>>("entropy_variable") = tmp;
    params.set<PostprocessorName>("pp_max_var") = _max_var_pp[i];
    params.set<PostprocessorName>("pp_min_var") = _min_var_pp[i];
    params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
    _problem->addAuxKernel(type, _aux_kernels[i], params);
  }
  for (unsigned i = 0; i < _temperature.size(); ++i)
  {
    tmp[0] = _temperature[i];
    params = _factory.getValidParams(type);
    _aux_kernels.push_back("entropy_" + _temperature[i] + "_aux");
    params.set<AuxVariableName>("variable") = _aux_variables.back();
    params.set<std::vector<VariableName>>("entropy_variable") = tmp;
    params.set<PostprocessorName>("pp_max_var") = _max_var_pp.back();
    params.set<PostprocessorName>("pp_min_var") = _min_var_pp.back();
    params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
    _problem->addAuxKernel(type, _aux_kernels.back(), params);
  }
}

/******************************************************************************/
/*                          Add Kernels                                       */
/******************************************************************************/
void
LynxAdvectionAction::createKernelActions()
{
  InputParameters params = emptyInputParameters();
  std::string name = "LynxAdvection";
  std::vector<VariableName> tmp;
  tmp.resize(1);
  for (unsigned i = 0; i < _compositional_phases.size(); ++i)
  {
    tmp[0] = _aux_variables[i];
    params = _factory.getValidParams(name + "Composition");
    params.set<NonlinearVariableName>("variable") = _compositional_phases[i];
    params.set<std::vector<VariableName>>("entropy") = tmp;
    params.set<std::vector<VariableName>>("velocities") = _velocities;
    params.set<std::vector<VariableName>>("displacements") = _displacements;
    params.set<MooseEnum>("element_length_type") = _element_length_type;
    params.set<Real>("beta_stabilization") = _beta_stabilization;
    params.set<Real>("cr_stabilization") = _cr_stabilization;
    params.set<PostprocessorName>("pp_max_vel") = "max_vel_pp";
    params.set<PostprocessorName>("pp_max_var") = _max_var_pp[i];
    params.set<PostprocessorName>("pp_min_var") = _min_var_pp[i];
    params.set<PostprocessorName>("pp_avg_var") = _avg_var_pp[i];
    params.set<PostprocessorName>("pp_max_entropy") = _max_entropy_pp[i];
    params.set<PostprocessorName>("pp_min_entropy") = _min_entropy_pp[i];
    params.set<PostprocessorName>("pp_avg_entropy") = _avg_entropy_pp[i];
    _problem->addKernel(name + "Composition", "advection_" + _compositional_phases[i], params);
    params = _factory.getValidParams("TimeDerivative");
    params.set<NonlinearVariableName>("variable") = _compositional_phases[i];
    _problem->addKernel("TimeDerivative", "time_" + _compositional_phases[i], params);
  }
  for (unsigned i = 0; i < _temperature.size(); ++i)
  {
    tmp[0] = _aux_variables.back();
    params = _factory.getValidParams(name + "Temperature");
    params.set<NonlinearVariableName>("variable") = _temperature[i];
    params.set<std::vector<VariableName>>("entropy") = tmp;
    params.set<std::vector<VariableName>>("velocities") = _velocities;
    params.set<std::vector<VariableName>>("displacements") = _displacements;
    params.set<MooseEnum>("element_length_type") = _element_length_type;
    params.set<Real>("beta_stabilization") = _beta_stabilization;
    params.set<Real>("cr_stabilization") = _cr_stabilization;
    params.set<Real>("coeff_shear_heating") = _coeff_Hs;
    params.set<PostprocessorName>("pp_max_vel") = "max_vel_pp";
    params.set<PostprocessorName>("pp_max_var") = _max_var_pp.back();
    params.set<PostprocessorName>("pp_min_var") = _min_var_pp.back();
    params.set<PostprocessorName>("pp_avg_var") = _avg_var_pp.back();
    params.set<PostprocessorName>("pp_max_entropy") = _max_entropy_pp.back();
    params.set<PostprocessorName>("pp_min_entropy") = _min_entropy_pp.back();
    params.set<PostprocessorName>("pp_avg_entropy") = _avg_entropy_pp.back();
    _problem->addKernel(name + "Temperature", "advection_" + _temperature[i], params);
    params = _factory.getValidParams("TimeDerivative");
    params.set<NonlinearVariableName>("variable") = _temperature[i];
    _problem->addKernel("TimeDerivative", "time_" + _temperature[i], params);
  }
}
