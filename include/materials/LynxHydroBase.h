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
