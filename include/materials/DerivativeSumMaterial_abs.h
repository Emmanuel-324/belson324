//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DerivativeFunctionMaterialBase.h"

template <bool is_ad>
class DerivativeSumMaterial_absTempl : public DerivativeFunctionMaterialBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  DerivativeSumMaterial_absTempl(const InputParameters & parameters);

  virtual void initialSetup();

protected:
  using DerivativeFunctionMaterialBaseTempl<is_ad>::_prop_F;
  using DerivativeFunctionMaterialBaseTempl<is_ad>::_qrule;
  using DerivativeFunctionMaterialBaseTempl<is_ad>::_qp;

  virtual void computeProperties();

  std::vector<std::string> _sum_materials;
  unsigned int _num_materials;

  /// User-defined ratio for gamma prime to gamma double prime
  Real _gamma_ratio;

  /// Prefactor for summation
  std::vector<Real> _prefactor;

  /// Function values of the summands.
  std::vector<const GenericMaterialProperty<Real, is_ad> *> _summand_F;
};

typedef DerivativeSumMaterial_absTempl<false> DerivativeSumMaterial_abs;
typedef DerivativeSumMaterial_absTempl<true> ADDerivativeSumMaterial_abs;
