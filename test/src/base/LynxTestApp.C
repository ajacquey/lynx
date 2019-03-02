//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "LynxTestApp.h"
#include "LynxApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

template <>
InputParameters
validParams<LynxTestApp>()
{
  InputParameters params = validParams<LynxApp>();
  return params;
}

LynxTestApp::LynxTestApp(InputParameters parameters) : MooseApp(parameters)
{
  LynxTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

LynxTestApp::~LynxTestApp() {}

void
LynxTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  LynxApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"LynxTestApp"});
    Registry::registerActionsTo(af, {"LynxTestApp"});
  }
}

void
LynxTestApp::registerApps()
{
  registerApp(LynxApp);
  registerApp(LynxTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
LynxTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  LynxTestApp::registerAll(f, af, s);
}
extern "C" void
LynxTestApp__registerApps()
{
  LynxTestApp::registerApps();
}
