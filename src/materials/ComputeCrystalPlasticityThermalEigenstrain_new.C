#include "ComputeCrystalPlasticityThermalEigenstrain_new.h"
#include "Function.h" // Include the full definition of Function

registerMooseObject("belson324App", ComputeCrystalPlasticityThermalEigenstrain_new);

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
        "_lattice_thermal_expansion_coefficients"))
{
  // Check if dynamic coefficients are used
  if (!_thermal_expansion_coeff_function_names.empty())
  {
    for (const std::string & name : _thermal_expansion_coeff_function_names)
    {
      _thermal_expansion_coeff_functions.push_back(&getFunction(name)); // Store pointers
    }
  }
}

void ComputeCrystalPlasticityThermalEigenstrain_new::computeQpDeformationGradient()
{
  RankTwoTensor thermal_coefficients = _thermal_expansion_coefficients;

  // If dynamic coefficients are provided, compute them
  if (!_thermal_expansion_coeff_function_names.empty())
  {
    for (size_t i = 0; i < _thermal_expansion_coeff_functions.size(); ++i)
    {
      Real value = _thermal_expansion_coeff_functions[i]->value(_t, _q_point[_qp]);
      thermal_coefficients(i, i) = value; // Update diagonal terms for anisotropic cases
    }
  }

  // Rotate the thermal deformation gradient for crystals based on Euler angles
  _lattice_thermal_expansion_coefficients[_qp] =
      _crysrot[_qp] * thermal_coefficients * _crysrot[_qp].transpose();

  // Compute the deformation gradient due to thermal expansion
  mooseAssert(_dt > 0, "Time increment (_dt) must be greater than zero.");
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
