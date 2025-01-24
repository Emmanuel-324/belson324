//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeCrystalPlasticityThermalEigenstrain_new.h"

registerMooseObject("belson324App", ComputeCrystalPlasticityThermalEigenstrain_new);

InputParameters
ComputeCrystalPlasticityThermalEigenstrain_new::validParams()
{
  InputParameters params = ComputeCrystalPlasticityEigenstrainBase::validParams();
  params.addClassDescription("Computes the deformation gradient associated with the linear thermal "
                             "expansion in a crystal plasticity simulation");
  params.addCoupledVar("temperature", "Coupled temperature variable");

  // Add new parameters for dynamic coefficients
  params.addParam<std::vector<std::string>>(
      "thermal_expansion_coeff_function_names",
      "Optional names of auxiliary variables or functions defining dynamic coefficients.");

  // Keep the existing fixed coefficients parameter
  params.addRequiredRangeCheckedParam<std::vector<Real>>(
      "thermal_expansion_coefficients",
      "thermal_expansion_coefficients_size=1 | thermal_expansion_coefficients_size=3 | "
      "thermal_expansion_coefficients_size=6 | thermal_expansion_coefficients_size=9",
      "Vector of values defining the constant second order thermal expansion coefficients, "
      "depending on the degree of anisotropy, this should be of size 1, 3, 6 or 9");

  return params;
}

ComputeCrystalPlasticityThermalEigenstrain_new::ComputeCrystalPlasticityThermalEigenstrain_new(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<ComputeCrystalPlasticityEigenstrainBase>(parameters),
    _temperature(coupledValue("temperature")),
    _temperature_old(coupledValueOld("temperature")),
    _ddeformation_gradient_dT(isCoupledConstant("temperature")
                                  ? nullptr
                                  : &declarePropertyDerivative<RankTwoTensor>(
                                        _deformation_gradient_name, coupledName("temperature", 0))),
    _thermal_expansion_coefficients(getParam<std::vector<Real>>("thermal_expansion_coefficients")),
    _thermal_expansion_coeff_function_names(getParam<std::vector<std::string>>(
        "thermal_expansion_coeff_function_names")),
    _lattice_thermal_expansion_coefficients(declareProperty<RankTwoTensor>(
        _eigenstrain_name +
        "_lattice_thermal_expansion_coefficients")) // avoid duplicated material name by including
                                                    // the eigenstrain name this coeff corresponds
                                                    // to
{
  // Check if dynamic coefficients are used
  if (!_thermal_expansion_coeff_function_names.empty())
  {
    for (const std::string & name : _thermal_expansion_coeff_function_names)
    {
      _thermal_expansion_coeff_functions.push_back(getFunction(name));
    }
  }
}


void
ComputeCrystalPlasticityThermalEigenstrain_new::computeQpDeformationGradient()
{
  std::vector<Real> thermal_coefficients = _thermal_expansion_coefficients;

  // If dynamic coefficients are provided, compute them
  if (!_thermal_expansion_coeff_function_names.empty())
  {
    for (size_t i = 0; i < _thermal_expansion_coeff_functions.size(); ++i)
    {
      thermal_coefficients[i] =
          _thermal_expansion_coeff_functions[i]->value(_t, _q_point[_qp]);
    }
  }

  // Rotate the thermal deformation gradient for crystals based on Euler angles
  _lattice_thermal_expansion_coefficients[_qp] =
      thermal_coefficients.rotated(_crysrot[_qp]);

  // Compute the deformation gradient due to thermal expansion
  Real dtheta = (_temperature[_qp] - _temperature_old[_qp]) * _substep_dt / _dt;
  RankTwoTensor residual_equivalent_thermal_expansion_increment =
      RankTwoTensor::Identity() - dtheta * _lattice_thermal_expansion_coefficients[_qp];
  _deformation_gradient[_qp] =
      residual_equivalent_thermal_expansion_increment.inverse() * _deformation_gradient_old[_qp];

  // Compute the derivative of deformation gradient w.r.t temperature
  if (_ddeformation_gradient_dT)
    (*_ddeformation_gradient_dT)[_qp] =
        _lattice_thermal_expansion_coefficients[_qp] * _deformation_gradient[_qp];
}
