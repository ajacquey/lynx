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

#include "LynxVelocityNormAux.h"

registerMooseObject("LynxApp", LynxVelocityNormAux);

template <>
InputParameters
validParams<LynxVelocityNormAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("velocities", "The string of velocities that advect the problem.");
  params.addClassDescription("Calculates the norm of the solid velocity vector.");
  return params;
}

LynxVelocityNormAux::LynxVelocityNormAux(const InputParameters & parameters)
  : AuxKernel(parameters), _n_vel(coupledComponents("velocities")), _vel(3)
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_n_vel != _mesh.dimension())
    mooseError("The number of variables supplied in 'velocities' must match the mesh dimension.");

  // Fetch coupled variables and gradients (as stateful properties if necessary)
  for (unsigned int i = 0; i < _n_vel; ++i)
    _vel[i] = &coupledValue("velocities", i);

  // Set unused dimensions to zero
  for (unsigned i = _n_vel; i < 3; ++i)
    _vel[i] = &_zero;
}

Real
LynxVelocityNormAux::computeValue()
{
  RealVectorValue v((*_vel[0])[_qp], (*_vel[1])[_qp], (*_vel[2])[_qp]);

  return v.norm();
}
