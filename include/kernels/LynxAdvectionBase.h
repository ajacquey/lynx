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

#ifndef LYNXADVECTIONBASE_H
#define LYNXADVECTIONBASE_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

class LynxAdvectionBase;

template <>
InputParameters validParams<LynxAdvectionBase>();

class LynxAdvectionBase : public DerivativeMaterialInterface<Kernel>
{
public:
  LynxAdvectionBase(const InputParameters & parameters);
  ~LynxAdvectionBase() {}
  // enum ElementLengthType
  // {
  //   MIN,
  //   MAX,
  //   AVERAGE
  // };

protected:
  virtual void precalculateResidual() override;
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int) override;

  virtual Real computeElementDiameter();
  virtual Real computeArtificialViscosity() = 0;

  MooseEnum _element_length_type;

  Real _beta_stabilization;
  Real _cr_stabilization;

  // velocities
  bool _has_disp;
  unsigned int _nvel;
  std::vector<unsigned> _disp_var;
  std::vector<const VariableValue *> _vel;
  std::vector<const VariableValue *> _vel_old;
  std::vector<const VariableValue *> _vel_older;

  //  old advected quantity
  const VariableValue & _value_old;
  const VariableGradient & _gradient_old;
  const VariableSecond & _second_old;

  //  older advected quantity
  const VariableValue & _value_older;
  const VariableGradient & _gradient_older;
  const VariableSecond & _second_older;

  // old and older entropy
  const VariableValue & _entropy_old;
  const VariableValue & _entropy_older;

  // postprocessor to get the global maxima and minima of all quantities
  const PostprocessorValue & _pp_max_vel;
  const PostprocessorValue & _pp_max_var;
  const PostprocessorValue & _pp_min_var;
  const PostprocessorValue & _pp_avg_var;
  const PostprocessorValue & _pp_max_entropy;
  const PostprocessorValue & _pp_min_entropy;
  const PostprocessorValue & _pp_avg_entropy;

  std::vector<Real> _residual;
  Real _artificial_viscosity;
};

#endif // LYNXADVECTIONBASE_H
