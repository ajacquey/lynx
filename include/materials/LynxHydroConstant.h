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
