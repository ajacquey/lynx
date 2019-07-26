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
