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

#ifndef LYNXSTRAINRATEAUX_H
#define LYNXSTRAINRATEAUX_H

#include "LynxStrainAuxBase.h"

class LynxStrainRateAux;

template <>
InputParameters validParams<LynxStrainAux>();

class LynxStrainRateAux : public LynxStrainAuxBase
{
public:
  LynxStrainRateAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:
  const unsigned int _i;
  const unsigned int _j;
};

#endif // LYNXSTRAINRATEAUX_H
