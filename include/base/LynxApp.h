//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef LYNXAPP_H
#define LYNXAPP_H

#include "MooseApp.h"

class LynxApp;

template <>
InputParameters validParams<LynxApp>();

class LynxApp : public MooseApp
{
public:
  LynxApp(InputParameters parameters);
  virtual ~LynxApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};

#endif /* LYNXAPP_H */
