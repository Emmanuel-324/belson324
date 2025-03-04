//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DerivativeParsedMaterial_abs.h"

registerMooseObject("belson324App", DerivativeParsedMaterial_abs);
registerMooseObject("belson324App", ADDerivativeParsedMaterial_abs);

template <bool is_ad>
InputParameters
DerivativeParsedMaterial_absTempl<is_ad>::validParams()
{
  InputParameters params = DerivativeParsedMaterialHelperTempl<is_ad>::validParams();
  params += ParsedMaterialBase::validParams();
  params.addClassDescription("Parsed Function Material with automatic derivatives.");
  
  // Adding gamma_ratio as a user-defined parameter
  params.addParam<Real>("gamma_ratio", 0.25, "Ratio of gamma prime to gamma double prime.");
  
  return params;
}

template <bool is_ad>
DerivativeParsedMaterial_absTempl<is_ad>::DerivativeParsedMaterial_absTempl(
    const InputParameters & parameters)
  : DerivativeParsedMaterialHelperTempl<is_ad>(parameters,
                                               Moose::VariableNameMappingMode::USE_MOOSE_NAMES),
    ParsedMaterialBase(parameters, this),
    _gamma_ratio(this->template getParam<Real>("gamma_ratio")) // Store user-defined ratio
{
  // Adjust free energy function to include user-defined gamma_ratio
  Real constrained_c1 = _gamma_ratio * this->coupledValue("eta2"); // Corrected eta2 reference

  // Build function, take derivatives, optimize
  this->functionParse(_function, _constant_names, _constant_expressions, constrained_c1,
                      _tol_names, _tol_values, _functor_names, _functor_symbols);
}

// explicit instantiation
template class DerivativeParsedMaterial_absTempl<false>;
template class DerivativeParsedMaterial_absTempl<true>;
