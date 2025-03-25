//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

// Forward Declarations

/**
 * Enforce the equality of the chemical potentials in a ternary system.
 * Eq. (21) in the original KKS paper extended to two solutes.
 *
 * \f$ dF_{\gamma}/dc_{\text{Al}} = dF_{\gamma'}/dc_{\text{Al}} \f$
 * \f$ dF_{\gamma}/dc_{\text{Nb}} = dF_{\gamma'}/dc_{\text{Nb}} \f$
 *
 * We supply free energy functions (i.e., KKSBaseMaterial) for each solute.
 *
 * \see KKSPhaseConcentrationTernary
 */
class KKSPhaseChemicalPotentialTernary
  : public DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>
{
public:
  static InputParameters validParams();

  KKSPhaseChemicalPotentialTernary(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  virtual void initialSetup();

private:
  /// Coupled variables for solutes
  unsigned int _c_Nb_var;
  VariableName _c_Nb_name;

  /// Material properties for free energy derivatives
  const MaterialProperty<Real> & _dfgammadc_Al;
  const MaterialProperty<Real> & _dfgammaPdc_Al;
  const MaterialProperty<Real> & _dfgammadc_Nb;
  const MaterialProperty<Real> & _dfgammaPdc_Nb;

  const MaterialProperty<Real> & _d2fgammadc_Al2;
  const MaterialProperty<Real> & _d2fgammaPdc_Al2;
  const MaterialProperty<Real> & _d2fgammadc_Nb2;
  const MaterialProperty<Real> & _d2fgammaPdc_Nb2;

  std::vector<const MaterialProperty<Real> *> _d2fgammadc_Al_darg;
  std::vector<const MaterialProperty<Real> *> _d2fgammadc_Nb_darg;

  /// Site fractions
  const Real _k_Al;
  const Real _k_Nb;
};
