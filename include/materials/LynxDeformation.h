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

#ifndef LYNXDEFORMATION_H
#define LYNXDEFORMATION_H

#include "LynxDeformationBase.h"

class LynxDeformation;

template <>
InputParameters validParams<LynxDeformation>();

class LynxDeformation : public LynxDeformationBase
{
public:
  LynxDeformation(const InputParameters & parameters);
  virtual ~LynxDeformation() {}

protected:
  virtual void initQpStatefulProperties() override;
  virtual void plasticCorrection(Real & pressure, RankTwoTensor & stress_dev) override;
  virtual Real computePlasticityYield(const Real & pressure, const Real & eqv_stress);
  virtual Real plasticIncrement(const Real & /*pressure*/, const Real & eqv_stress);
  virtual void computePlasticityProperties(const Real & pressure);
  virtual void updatePlasticityParameters();

  // Plastic parameters
  const std::vector<Real> _friction_angle_0;
  const std::vector<Real> _cohesion_0;
  const std::vector<Real> _friction_angle_res;
  const std::vector<Real> _cohesion_res;
  const std::vector<Real> _dilation_angle;
  const std::vector<Real> _intnl_0;
  const std::vector<Real> _intnl_lim;
  std::vector<Real> _one_on_plastic_eta;
  const bool _has_hardening;

  // Plasticity structure
  plasticity * _plasticity;

  // Plastic properties
  MaterialProperty<Real> * _intnl;
  const MaterialProperty<Real> * _intnl_old;
};

#endif // LYNXDEFORMATION_H
