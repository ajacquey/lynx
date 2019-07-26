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

#ifndef LYNXHYDROBASE_H
#define LYNXHYDROBASE_H

#include "LynxMaterialBase.h"
#include "RankTwoTensor.h"

class LynxHydroBase;

template <>
InputParameters validParams<LynxHydroBase>();

class LynxHydroBase : public LynxMaterialBase
{
public:
  LynxHydroBase(const InputParameters & parameters);
  virtual ~LynxHydroBase() {}

protected:
  virtual void computeQpProperties() override;
  virtual void computeQpCompressibilities();
  virtual void computeQpFluidMobility();
  virtual void computeQpPoroMech();
  virtual void computeQpFluidCompressibility() = 0;
  virtual void computeQpSolidCompressibility() = 0;
  virtual void computeQpPermeability() = 0;
  virtual void computeQpFluidViscosity() = 0;

  const VariableValue & _porosity;
  const MaterialProperty<Real> & _K;
  const MaterialProperty<RankFourTensor> & _tangent_modulus;
  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _viscous_strain_incr;
  const MaterialProperty<RankTwoTensor> & _plastic_strain_incr;
  MaterialProperty<Real> & _biot;
  MaterialProperty<Real> & _C_d;
  MaterialProperty<Real> & _C_biot;
  MaterialProperty<Real> & _fluid_mobility;
  MaterialProperty<Real> & _poro_mech;
  MaterialProperty<Real> & _poro_mech_jac;

  std::vector<Real> _C_f;
  std::vector<Real> _C_s;
  std::vector<Real> _k;
  std::vector<Real> _eta_f;
};

#endif // LYNXHYDROBASE_H
