#pragma once

#include "KernelValue.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

// Forward Declarations

/**
 * Enforce sum of phase concentrations to be the real concentration for a ternary system.
 *
 * \f$ c_{Al} = h_1 c_{Al}^{(1)} + h_2 c_{Al}^{(2)} + h_3 c_{Al}^{(3)} \f$
 * \f$ c_{Nb} = h_1 c_{Nb}^{(1)} + h_2 c_{Nb}^{(2)} + h_3 c_{Nb}^{(3)} \f$
 *
 * The non-linear variable for this Kernel is one of  the concentrations \f$ c_i \f$, while
 * \f$ c_j \neq c_i \f$ and \f$ c \f$ are supplied as coupled variables.
 *
 * \see KKSPhaseChemicalPotential
 */
class KKSMultiPhaseConcentrationTernary
  : public DerivativeMaterialInterface<JvarMapKernelInterface<KernelValue>>
{
public:
  static InputParameters validParams();

  KKSMultiPhaseConcentrationTernary(const InputParameters & parameters);

protected:
  virtual Real precomputeQpResidual();
  virtual Real precomputeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  // Number of phases
  const unsigned int _num_j;

  // Phase concentrations for each solute (Al, Nb)
  const std::vector<const VariableValue *> _cj_Al;
  const std::vector<const VariableValue *> _cj_Nb;

  // Mapping for coupled variables
  const JvarMap & _cj_Al_map;
  const JvarMap & _cj_Nb_map;

  // Position of the nonlinear variable in the list of cj's
  int _k;

  // Total concentrations of Al and Nb
  const VariableValue & _c_Al;
  const VariableValue & _c_Nb;

  // Switching functions for Al and Nb
  std::vector<MaterialPropertyName> _hj_Al_names;
  std::vector<MaterialPropertyName> _hj_Nb_names;

  // Material properties for interpolation functions
  std::vector<const MaterialProperty<Real> *> _prop_hj_Al;
  std::vector<const MaterialProperty<Real> *> _prop_hj_Nb;

  // Order parameters for each phase
  std::vector<VariableName> _eta_names;
  const JvarMap & _eta_map;

  // Derivative of the switching function wrt order parameters
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_dhjdetai;
};
