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

#ifndef LYNXHYDROCONSTANT
#define LYNXHYDROCONSTANT

#include "LynxHydroBase.h"

class LynxHydroConstant;

template <>
InputParameters validParams<LynxHydroConstant>();

class LynxHydroConstant : public LynxHydroBase
{
public:
  LynxHydroConstant(const InputParameters & parameters);

protected:
  virtual void computeQpFluidCompressibility() override;
  virtual void computeQpSolidCompressibility() override;
  virtual void computeQpPermeability() override;
  virtual void computeQpFluidViscosity() override;

  std::vector<Real> _perm;
  std::vector<Real> _fluid_viscosity;
  std::vector<Real> _fluid_compr;
  std::vector<Real> _solid_compr;
};

#endif // LYNXHYDROCONSTANT