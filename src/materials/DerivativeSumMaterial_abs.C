//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DerivativeSumMaterial_abs.h"
#include "libmesh/quadrature.h"

registerMooseObject("belson324App", DerivativeSumMaterial_abs);
registerMooseObject("belson324App", ADDerivativeSumMaterial_abs);

template <bool is_ad>
InputParameters
DerivativeSumMaterial_absTempl<is_ad>::validParams()
{
  InputParameters params = DerivativeFunctionMaterialBaseTempl<is_ad>::validParams();
  params.addClassDescription("Meta-material to sum up multiple derivative materials");
  params.addParam<std::vector<std::string>>("sum_materials",
                                            "Base name of the parsed sum material property");

  // User-defined gamma ratio
  params.addParam<Real>("gamma_ratio", 0.25, "Ratio of gamma prime to gamma double prime.");

  params.addRequiredCoupledVar("args", "Vector of names of variables being summed");
  params.addCoupledVar("displacement_gradients",
                       "Vector of displacement gradient variables (see "
                       "Modules/PhaseField/DisplacementGradients "
                       "action)");

  return params;
}

template <bool is_ad>
DerivativeSumMaterial_absTempl<is_ad>::DerivativeSumMaterial_absTempl(const InputParameters & parameters)
  : DerivativeFunctionMaterialBaseTempl<is_ad>(parameters),
    _sum_materials(this->template getParam<std::vector<std::string>>("sum_materials")),
    _num_materials(_sum_materials.size()),
    _gamma_ratio(this->template getParam<Real>("gamma_ratio")),
    _prefactor(_num_materials, 1.0)
{
  if (_num_materials == 0)
    mooseError("Please supply at least one material to sum in DerivativeSumMaterial_abs", name());

  for (unsigned int n = 0; n < _num_materials; ++n)
  {
    _summand_F.push_back(&this->template getGenericMaterialProperty<Real, is_ad>(_sum_materials[n]));
  }
}

template <bool is_ad>
void
DerivativeSumMaterial_absTempl<is_ad>::computeProperties()
{
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    (*_prop_F)[_qp] = (*_summand_F[0])[_qp] * _gamma_ratio; // Apply user-defined ratio
    for (unsigned int n = 1; n < _num_materials; ++n)
    {
      (*_prop_F)[_qp] += (*_summand_F[n])[_qp];
    }
  }
}

// Explicit instantiation
template class DerivativeSumMaterial_absTempl<false>;
template class DerivativeSumMaterial_absTempl<true>;
