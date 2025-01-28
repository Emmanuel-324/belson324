#pragma once

#include "ComputeCrystalPlasticityEigenstrainBase.h"
#include "DerivativeMaterialInterface.h"

/**
 * ComputeCrystalPlasticityThermalEigenstrain_new computes an eigenstrain for thermal expansion
 * with a constant or variable thermal expansion coefficient.
 * The formulations for thermal deformation gradient and thermal eigenstrain are as follows:
 *    \dot{Fth} * Fth^{-1} = \dot{_temperature} * _thermal_expansion_coefficients;
 *    thermal_eigenstrain = 0.5 * (Fth^{T} * Fth - I)
 */
class ComputeCrystalPlasticityThermalEigenstrain_new
  : public DerivativeMaterialInterface<ComputeCrystalPlasticityEigenstrainBase>
{
public:
  static InputParameters validParams();

  ComputeCrystalPlasticityThermalEigenstrain_new(const InputParameters & parameters);

protected:
  /// Compute the deformation gradient due to thermal expansion
  virtual void computeQpDeformationGradient() override;

  /// Temperature variable value
  const VariableValue & _temperature;
  const VariableValue & _temperature_old;

  /// Stores the derivative of the deformation gradient w.r.t. temperature
  MaterialProperty<RankTwoTensor> * _ddeformation_gradient_dT;

  /// The thermal expansion coefficient defined in the Cartesian coordinate system
  const RankTwoTensor _thermal_expansion_coefficients;

  /// Stores the thermal expansion coefficient w.r.t. the lattice symmetry axis
  MaterialProperty<RankTwoTensor> & _lattice_thermal_expansion_coefficients;

  /// Names of the thermal expansion coefficient functions (if variable)
  std::vector<std::string> _thermal_expansion_coeff_function_names;

  /// Functions for variable thermal expansion coefficients
  std::vector<Function *> _thermal_expansion_coeff_functions;
};
