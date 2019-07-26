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

#ifndef LYNXVOLSTRAINRATEAUX_H
#define LYNXVOLSTRAINRATEAUX_H

#include "LynxStrainAuxBase.h"

class LynxVolStrainRateAux;

template <>
InputParameters validParams<LynxVolStrainRateAux>();

class LynxVolStrainRateAux : public LynxStrainAuxBase
{
public:
  LynxVolStrainRateAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
};

#endif // LYNXVOLSTRAINRATEAUX_H
