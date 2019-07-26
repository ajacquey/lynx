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

#ifndef LYNXSTRAINAUX_H
#define LYNXSTRAINAUX_H

#include "LynxStrainAuxBase.h"

class LynxStrainAux;

template <>
InputParameters validParams<LynxStrainAux>();

class LynxStrainAux : public LynxStrainAuxBase
{
public:
  LynxStrainAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:
  const unsigned int _i;
  const unsigned int _j;
};

#endif // LYNXSTRAINAUX_H
