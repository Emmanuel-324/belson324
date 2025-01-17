//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BarrierFunctionMaterial_new.h"

registerMooseObject("belson324App", BarrierFunctionMaterial_new);

InputParameters
BarrierFunctionMaterial_new::validParams()
{
  InputParameters params = OrderParameterFunctionMaterial::validParams();
  params.addClassDescription("Helper material to provide $g(\\eta)$ and its derivative in a "
                             "polynomial.\nSIMPLE: $\\eta^2(1-\\eta)^2$\nLOW: $\\eta(1-\\eta)$"
                             "\nHIGH: $\\eta^2(1-\\eta^2)^2$");
  MooseEnum g_order("SIMPLE=0 LOW HIGH", "SIMPLE");
  params.addParam<MooseEnum>("g_order", g_order, "Polynomial order of the barrier function g(eta)");
  params.addParam<bool>("well_only",
                        false,
                        "Make the g zero in [0:1] so it only contributes to "
                        "enforcing the eta range and not to the phase "
                        "transformation barrier.");
  params.set<std::string>("function_name") = std::string("g");
  return params;
}

BarrierFunctionMaterial_new::BarrierFunctionMaterial_new(const InputParameters & parameters)
  : OrderParameterFunctionMaterial(parameters),
    _g_order(getParam<MooseEnum>("g_order")),
    _well_only(getParam<bool>("well_only"))
{
}

void
BarrierFunctionMaterial_new::computeQpProperties()
{
  // Clamping eta to the range [0, 1]
  Real clamped_eta = std::max(0.0, std::min(_eta[_qp], 1.0));

  if (_well_only && clamped_eta >= 0.0 && clamped_eta <= 1.0)
  {
    _prop_f[_qp] = 0.0;
    _prop_df[_qp] = 0.0;
    _prop_d2f[_qp] = 0.0;
    return;
  }

  switch (_g_order)
  {
    case 0: // SIMPLE
      _prop_f[_qp] = clamped_eta * clamped_eta * (1.0 - clamped_eta) * (1.0 - clamped_eta);
      _prop_df[_qp] = 2.0 * clamped_eta * (clamped_eta - 1.0) * (2.0 * clamped_eta - 1.0);
      _prop_d2f[_qp] = 12.0 * (clamped_eta * clamped_eta - clamped_eta) + 2.0;
      break;

    case 1: // LOW
      _prop_f[_qp] = clamped_eta * (1.0 - clamped_eta);
      _prop_df[_qp] = 1.0 - 2.0 * clamped_eta;
      _prop_d2f[_qp] = -2.0;
      break;

    case 2: // HIGH
      _prop_f[_qp] = clamped_eta * clamped_eta * (1.0 - clamped_eta * clamped_eta) * (1.0 - clamped_eta * clamped_eta);
      _prop_df[_qp] = clamped_eta * (2.0 - clamped_eta * clamped_eta * (8.0 + 6.0 * clamped_eta * clamped_eta));
      _prop_d2f[_qp] = 2.0 - clamped_eta * clamped_eta * (24.0 + 30.0 * clamped_eta * clamped_eta);
      break;

    default:
      mooseError("Internal error");
  }
}
