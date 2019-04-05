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

static void
associateSyntaxInner(Syntax & syntax, ActionFactory & /*action_factory*/)
{
  registerSyntax("LynxAdvectionAction", "LynxAdvection");

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
