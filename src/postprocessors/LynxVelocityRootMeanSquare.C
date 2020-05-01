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

#include "LynxVelocityRootMeanSquare.h"
#include "Function.h"

registerMooseObject("MooseApp", LynxVelocityRootMeanSquare);

InputParameters
LynxVelocityRootMeanSquare::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Calculates the mean square root of the velocity vector.");
  params.addRequiredCoupledVar("velocities", "The string of velocities that advect the problem.");
  return params;
}

LynxVelocityRootMeanSquare::LynxVelocityRootMeanSquare(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters), _n_vel(coupledComponents("velocities")), _vel(3)
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
LynxVelocityRootMeanSquare::getValue()
{
  return std::sqrt(ElementIntegralPostprocessor::getValue());
}

Real
LynxVelocityRootMeanSquare::computeQpIntegral()
{
  RealVectorValue vel((*_vel[0])[_qp], (*_vel[1])[_qp], (*_vel[2])[_qp]);

  return vel.norm_sq();
}
