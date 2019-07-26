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
