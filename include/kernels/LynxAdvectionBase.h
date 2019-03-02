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
