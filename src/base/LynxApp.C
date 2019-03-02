#include "LynxApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<LynxApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

LynxApp::LynxApp(InputParameters parameters) : MooseApp(parameters)
{
  LynxApp::registerAll(_factory, _action_factory, _syntax);
}

LynxApp::~LynxApp() {}

void
LynxApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"LynxApp"});
  Registry::registerActionsTo(af, {"LynxApp"});

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
