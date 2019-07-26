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

#ifndef LYNXADVECTIONCOMPOSITION_H
#define LYNXADVECTIONCOMPOSITION_H

#include "LynxAdvectionBase.h"

class LynxAdvectionComposition;

template <>
InputParameters validParams<LynxAdvectionComposition>();

class LynxAdvectionComposition : public LynxAdvectionBase
{
public:
  LynxAdvectionComposition(const InputParameters & parameters);

protected:
  virtual Real computeArtificialViscosity() override;
  virtual void computeEntropyResidual();
};

#endif // LYNXADVECTIONCOMPOSITION_H
