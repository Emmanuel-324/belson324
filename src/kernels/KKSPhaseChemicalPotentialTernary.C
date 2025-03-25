//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSPhaseChemicalPotentialTernary.h"
#include "MathUtils.h"

using namespace MathUtils;

registerMooseObject("belson324App", KKSPhaseChemicalPotentialTernary);

InputParameters
KKSPhaseChemicalPotentialTernary::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Ternary KKS kernel to enforce the equality of chemical potentials "
                             "between phases for multiple solutes.");

  params.addRequiredCoupledVar("c_Nb", "Concentration of Nb"); 
  params.addRequiredParam<MaterialPropertyName>("fgamma_name",
                                                "Base name of the free energy function "
                                                "F_gamma");
  params.addRequiredParam<MaterialPropertyName>("fgammaP_name",
                                                "Base name of the free energy function "
                                                "F_gamma'");

  params.addParam<Real>("k_Al",
                        1.0,
                        "Site fraction for Al (e.g., sublattice concentration)");

  params.addParam<Real>("k_Nb",
                        1.0,
                        "Site fraction for Nb (e.g., sublattice concentration)");

  return params;
}

KKSPhaseChemicalPotentialTernary::KKSPhaseChemicalPotentialTernary(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    _c_Nb_var(coupled("c_Nb")),
    _c_Nb_name(coupledName("c_Nb", 0)),
    // First derivatives
    _dfgammadc_Al(getMaterialPropertyDerivative<Real>("fgamma_name", _var.name())),
    _dfgammaPdc_Al(getMaterialPropertyDerivative<Real>("fgammaP_name", _var.name())),
    _dfgammadc_Nb(getMaterialPropertyDerivative<Real>("fgamma_name", _c_Nb_name)),
    _dfgammaPdc_Nb(getMaterialPropertyDerivative<Real>("fgammaP_name", _c_Nb_name)),
    // Second derivatives
    _d2fgammadc_Al2(getMaterialPropertyDerivative<Real>("fgamma_name", _var.name(), _var.name())),
    _d2fgammaPdc_Al2(getMaterialPropertyDerivative<Real>("fgammaP_name", _var.name(), _var.name())),
    _d2fgammadc_Nb2(getMaterialPropertyDerivative<Real>("fgamma_name", _c_Nb_name, _c_Nb_name)),
    _d2fgammaPdc_Nb2(getMaterialPropertyDerivative<Real>("fgammaP_name", _c_Nb_name, _c_Nb_name)),
    _d2fgammadc_Al_darg(_n_args),
    _d2fgammadc_Nb_darg(_n_args),
    // Site fractions
    _k_Al(getParam<Real>("k_Al")),
    _k_Nb(getParam<Real>("k_Nb"))
{
  for (std::size_t i = 0; i < _n_args; ++i)
  {
    _d2fgammadc_Al_darg[i] = &getMaterialPropertyDerivative<Real>("fgamma_name", _var.name(), i);
    _d2fgammadc_Nb_darg[i] = &getMaterialPropertyDerivative<Real>("fgamma_name", _c_Nb_name, i);
  }
}

void
KKSPhaseChemicalPotentialTernary::initialSetup()
{
  validateNonlinearCoupling<Real>("fgamma_name");
  validateNonlinearCoupling<Real>("fgammaP_name");
}

Real
KKSPhaseChemicalPotentialTernary::computeQpResidual()
{
  return _test[_i][_qp] * (_dfgammadc_Al[_qp] / _k_Al - _dfgammaPdc_Al[_qp] / _k_Al +
                           _dfgammadc_Nb[_qp] / _k_Nb - _dfgammaPdc_Nb[_qp] / _k_Nb);
}

Real
KKSPhaseChemicalPotentialTernary::computeQpJacobian()
{
  return _test[_i][_qp] * _phi[_j][_qp] *
         (_d2fgammadc_Al2[_qp] / _k_Al - _d2fgammaPdc_Al2[_qp] / _k_Al +
          _d2fgammadc_Nb2[_qp] / _k_Nb - _d2fgammaPdc_Nb2[_qp] / _k_Nb);
}

Real
KKSPhaseChemicalPotentialTernary::computeQpOffDiagJacobian(unsigned int jvar)
{
  const unsigned int cvar = mapJvarToCvar(jvar);
  return _test[_i][_qp] * _phi[_j][_qp] *
         ((*_d2fgammadc_Al_darg[cvar])[_qp] / _k_Al - (*_d2fgammadc_Nb_darg[cvar])[_qp] / _k_Nb);
}
