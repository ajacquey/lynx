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


#include "LynxHoldStressAction.h"
#include "Factory.h"
#include "FEProblem.h"
#include "Conversion.h"

registerMooseAction("LynxApp", LynxHoldStressAction, "add_bc");

InputParameters
LynxHoldStressAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Set up hold stress boundary conditions.");
  params.addRequiredParam<std::vector<BoundaryName>>(
      "boundary", "The list of boundary IDs from the mesh where the pressure will be applied.");
  params.addRequiredParam<std::vector<NonlinearVariableName>>(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  // Elastic moduli parameters
  params.addRangeCheckedParam<Real>(
      "bulk_modulus", "bulk_modulus >= 0.0", "The drained bulk modulus of the material.");
  params.addRangeCheckedParam<Real>(
      "shear_modulus", "shear_modulus >= 0.0", "The shear modulus of the material.");
  params.addCoupledVar("fluid_pressure", "The fluid pressure variable.");
  return params;
}

LynxHoldStressAction::LynxHoldStressAction(const InputParameters & params) : Action(params) {}

void
LynxHoldStressAction::act()
{
  const std::string kernel_name = "LynxHoldStressBC";

  std::vector<NonlinearVariableName> displacements =
      getParam<std::vector<NonlinearVariableName>>("displacements");

  // Create pressure BCs
  for (unsigned int i = 0; i < displacements.size(); ++i)
  {
    // Create unique kernel name for each of the components
    std::string unique_kernel_name = kernel_name + "_" + _name + "_" + Moose::stringify(i);

    InputParameters params = _factory.getValidParams(kernel_name);
    params.applyParameters(parameters());
    params.set<bool>("use_displaced_mesh") = true;
    params.set<unsigned int>("component") = i;
    params.set<NonlinearVariableName>("variable") = displacements[i];

    _problem->addBoundaryCondition(kernel_name, unique_kernel_name, params);
  }
}
