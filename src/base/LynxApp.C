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

#include "LynxApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
LynxApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  params.set<bool>("use_legacy_material_output") = false;
  
  return params;
}

LynxApp::LynxApp(InputParameters parameters) : MooseApp(parameters)
{
  LynxApp::registerAll(_factory, _action_factory, _syntax);
}

LynxApp::~LynxApp() {}

static void
associateSyntaxInner(Syntax & syntax, ActionFactory & /*action_factory*/)
{
  registerSyntax("LynxAdvectionAction", "LynxAdvection");
  registerSyntax("LynxADAdvectionAction", "LynxADAdvection");

  registerSyntax("EmptyAction", "BCs/LynxPressure");
  registerSyntax("LynxPressureAction", "BCs/LynxPressure/*");

  registerSyntax("EmptyAction", "BCs/LynxHoldStress");
  registerSyntax("LynxHoldStressAction", "BCs/LynxHoldStress/*");
}

void
LynxApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"LynxApp"});
  Registry::registerActionsTo(af, {"LynxApp"});
  associateSyntaxInner(s, af);
  /* register custom execute flags, action syntax, etc. here */
}

void
LynxApp::registerApps()
{
  registerApp(LynxApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
LynxApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  LynxApp::registerAll(f, af, s);
}
extern "C" void
LynxApp__registerApps()
{
  LynxApp::registerApps();
}
