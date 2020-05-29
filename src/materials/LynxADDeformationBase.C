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

#include "LynxADDeformationBase.h"
#include "Assembly.h"

InputParameters
LynxADDeformationBase::validParams()
{
  InputParameters params = LynxADMaterialBase::validParams();
  params.addClassDescription("Base class calculating the deformation of a material.");
  // Coupled variables
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system.");
  params.addCoupledVar("lithostatic_pressure", "The lithostatic pressure variable.");
  params.addCoupledVar("temperature", "The temperature variable.");
  params.addCoupledVar("temperature_dot",
                       "The time derivative of the temperature auxialiary variable.");
  // Strain parameters
  MooseEnum strain_model("small=0 finite=1", "small");
  params.addParam<MooseEnum>(
      "strain_model", strain_model, "The model to use to calculate the strain rate tensor.");
  params.addParam<bool>(
      "volumetric_locking_correction", false, "Flag to correct volumetric locking");
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

LynxADDeformationBase::LynxADDeformationBase(const InputParameters & parameters)
  : LynxADMaterialBase(parameters),
    // Coupled variables
    _ndisp(coupledComponents("displacements")),
    _grad_disp(3),
    _grad_disp_old(3),
    _plith(isCoupled("lithostatic_pressure") ? adCoupledValue("lithostatic_pressure")
                                             : adZeroValue()),
    _coupled_temp(isCoupled("temperature")),
    _temp_dot(_coupled_temp ? adCoupledDot("temperature") : adZeroValue()),
    _coupled_temp_aux(isCoupled("temperature_dot")),
    _temp_dot_aux(_coupled_temp_aux ? adCoupledValue("temperature_dot") : adZeroValue()),
    // Strain parameters
    _strain_model(getParam<MooseEnum>("strain_model")),
    _vol_locking_correction(getParam<bool>("volumetric_locking_correction")),
    _current_elem_volume(_assembly.elemVolume()),
    // Strain properties
    _strain_increment(declareADProperty<RankTwoTensor>("strain_increment")),
    _spin_increment(declareADProperty<RankTwoTensor>("spin_increment")),
    _thermal_exp(_coupled_temp || _coupled_temp_aux
                     ? &getMaterialProperty<Real>("thermal_expansion_coefficient")
                     : nullptr),
    // Stress properties
    _stress(declareADProperty<RankTwoTensor>("stress")),
    _inelastic_heat(declareADProperty<Real>("inelastic_heat"))
{
  if (getParam<bool>("use_displaced_mesh"))
    paramError("use_displaced_mesh",
               "The strain and stress calculator needs to run on the undisplaced mesh.");

  if (_coupled_temp && _coupled_temp_aux)
    mooseWarning("LynxADDeformationBase: you provided both 'temperature' and 'temperature_dot'!");
}

void
LynxADDeformationBase::initialSetup()
{
  displacementIntegrityCheck();
  // fetch coupled variables and gradients (as stateful properties if necessary)
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _grad_disp[i] = &adCoupledGradient("displacements", i);
    if (_fe_problem.isTransient())
      _grad_disp_old[i] = &coupledGradientOld("displacements", i);
    else
      _grad_disp_old[i] = &_grad_zero;
  }

  // set unused dimensions to zero
  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _grad_disp[i] = &adZeroGradient();
    _grad_disp_old[i] = &_grad_zero;
  }
}

void
LynxADDeformationBase::displacementIntegrityCheck()
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    paramError(
        "displacements",
        "The number of variables supplied in 'displacements' must match the mesh dimension.");
}

void
LynxADDeformationBase::computeProperties()
{
  computeStrainIncrement();
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
    computeQpProperties();
}

void
LynxADDeformationBase::computeQpProperties()
{
  computeQpDeformation();
  computeQpThermalSources();
}

void
LynxADDeformationBase::computeStrainIncrement()
{
  ADReal vol_strain_incr = 0.0;
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    ADRankTwoTensor grad_tensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);
    RankTwoTensor grad_tensor_old(
        (*_grad_disp_old[0])[_qp], (*_grad_disp_old[1])[_qp], (*_grad_disp_old[2])[_qp]);

    switch (_strain_model)
    {
      case 0: // SMALL STRAIN
        calculateSmallStrain(grad_tensor, grad_tensor_old);
        break;
      case 1: // FINITE STRAIN
        calculateFiniteStrain(grad_tensor, grad_tensor_old);
        break;
      default:
        mooseError("Unknown strain model. Specify 'small' or 'finite'!");
    }

    // Thermal strain correction
    if (_coupled_temp && !_coupled_temp_aux) 
      _strain_increment[_qp].addIa(-(*_thermal_exp)[_qp] * _temp_dot[_qp] * _dt / 3.0);
    else if (_coupled_temp_aux && !_coupled_temp)
      _strain_increment[_qp].addIa(-(*_thermal_exp)[_qp] * _temp_dot_aux[_qp] * _dt / 3.0);

    if (_vol_locking_correction)
      vol_strain_incr += _strain_increment[_qp].trace() * _JxW[_qp] * _coord[_qp];
  }

  if (_vol_locking_correction)
  {
    vol_strain_incr /= _current_elem_volume;

    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
    {
      ADReal trace = _strain_increment[_qp].trace();
      _strain_increment[_qp](0, 0) += (vol_strain_incr - trace) / 3.0;
      _strain_increment[_qp](1, 1) += (vol_strain_incr - trace) / 3.0;
      _strain_increment[_qp](2, 2) += (vol_strain_incr - trace) / 3.0;
    }
  }
}

void
LynxADDeformationBase::calculateSmallStrain(const ADRankTwoTensor & grad_tensor,
                                            const RankTwoTensor & grad_tensor_old)
{
  ADRankTwoTensor A = grad_tensor - grad_tensor_old;

  _strain_increment[_qp] = 0.5 * (A + A.transpose());
  _spin_increment[_qp] = 0.5 * (A - A.transpose());
}

void
LynxADDeformationBase::calculateFiniteStrain(const ADRankTwoTensor & grad_tensor,
                                             const RankTwoTensor & grad_tensor_old)
{
  ADRankTwoTensor F = grad_tensor;
  RankTwoTensor F_old = grad_tensor_old;
  F.addIa(1.0);
  F_old.addIa(1.0);

  // Increment gradient
  ADRankTwoTensor L = -F_old * F.inverse();
  L.addIa(1.0);

  _strain_increment[_qp] = 0.5 * (L + L.transpose());
  _spin_increment[_qp] = 0.5 * (L - L.transpose());
}

void
LynxADDeformationBase::computeQpDeformation()
{
  // Initialize deformation
  initializeQpDeformation();

  // Compute Stress tensor
  computeQpStress();
}

void
LynxADDeformationBase::initializeQpDeformation()
{
}

void
LynxADDeformationBase::computeQpStress()
{
  // Update the volumetric part of the deformation
  ADReal pressure = volumetricDeformation();

  // Update the deviatoric part of the deformation
  ADRankTwoTensor stress_dev = deviatoricDeformation(pressure);

  // Form the total stress tensor
  reformStressTensor(pressure, stress_dev);
}

void
LynxADDeformationBase::reformStressTensor(const ADReal & pressure,
                                          const ADRankTwoTensor & stress_dev)
{
  _stress[_qp] = stress_dev;
  _stress[_qp].addIa(-pressure);
}

ADRankTwoTensor
LynxADDeformationBase::spinRotation(const ADRankTwoTensor & tensor)
{
  return tensor + _spin_increment[_qp] * tensor.deviatoric() - tensor.deviatoric() * _spin_increment[_qp];
}
