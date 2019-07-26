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

#include "LynxVelocityAux.h"

registerMooseObject("LynxApp", LynxVelocityAux);

template <>
InputParameters
validParams<LynxVelocityAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("displacement",
                               "The displacement variable to calculate the velocity.");
  params.addClassDescription("Calculates a component of the solid velocity based on displacement.");
  return params;
}

LynxVelocityAux::LynxVelocityAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _disp(coupledValue("displacement")),
    _disp_old(_is_transient ? coupledValueOld("displacement") : _zero)
{
}

Real
LynxVelocityAux::computeValue()
{
  Real inv_dt = 1.0;
  if (_is_transient)
    inv_dt /= _dt;

  return (_disp[_qp] - _disp_old[_qp]) * inv_dt;
}
