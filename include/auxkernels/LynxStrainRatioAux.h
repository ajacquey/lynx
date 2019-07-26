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

#ifndef LYNXSTRAINRATIOAUX_H
#define LYNXSTRAINRATIOAUX_H

#include "LynxElasticStrainAuxBase.h"

class LynxStrainRatioAux;

template <>
InputParameters validParams<LynxStrainRatioAux>();

class LynxStrainRatioAux : public LynxElasticStrainAuxBase
{
public:
  LynxStrainRatioAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();
};

#endif // LYNXSTRAINRATIOAUX_H
