//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BarrierFunctionMaterial_abs.h"

registerMooseObject("belson324App", BarrierFunctionMaterial_abs);

InputParameters
BarrierFunctionMaterial_abs::validParams()
{
  InputParameters params = OrderParameterFunctionMaterial::validParams();
  params.addClassDescription("Helper material to provide g(eta) and its derivative in a "
                             "polynomial.\nSIMPLE: eta^2*(1-eta)^2\nLOW: eta*(1-eta)"
                             "\nHIGH: eta^2*(1-eta^2)^2");
  MooseEnum g_order("SIMPLE=0 LOW HIGH", "SIMPLE");
  params.addParam<MooseEnum>("g_order", g_order, "Polynomial order of the barrier function g(eta)");
  params.addParam<bool>("well_only",
                        false,
                        "Make the g zero in [0:1] so it only contributes to "
                        "enforcing the eta range and not to the phase "
                        "transformation berrier.");
  params.set<std::string>("function_name") = std::string("g");
  return params;
}

BarrierFunctionMaterial_abs::BarrierFunctionMaterial_abs(const InputParameters & parameters)
  : OrderParameterFunctionMaterial(parameters),
    _g_order(getParam<MooseEnum>("g_order")),
    _well_only(getParam<bool>("well_only"))
{
}

void
BarrierFunctionMaterial_abs::computeQpProperties()
{
  const Real n = std::abs(_eta[_qp]); //orginal _eta[_qp]

  int sign_eta = 0;

  if (_eta[_qp] > 0) sign_eta = 1;
  else if (_eta[_qp] < 0) sign_eta = -1;

  if (_well_only && n >= 0.0 && n <= 1.0)
  {
    _prop_f[_qp] = 0.0;
    _prop_df[_qp] = 0.0;
    _prop_d2f[_qp] = 0.0;
    return;
  }

  switch (_g_order)
  {
    case 0: // SIMPLE
      _prop_f[_qp] = n * n * (1.0 - n) * (1.0 - n);
      _prop_df[_qp] = 2.0 * n * (n - 1.0) * (2.0 * n - 1.0) * sign_eta;
      _prop_d2f[_qp] = (12.0 * (n * n - n) + 2.0) * sign_eta;
      break;

    case 1: // LOW
      _prop_f[_qp] = n * (1.0 - n);
      _prop_df[_qp] = 1.0 - 2.0 * n * sign_eta;
      _prop_d2f[_qp] = -2.0 * sign_eta;
      break;

    case 2: // HIGH
      _prop_f[_qp] = n * n * (1.0 - n * n) * (1.0 - n * n);
      _prop_df[_qp] = n * (2.0 - n * n * (8.0 + 6.0 * n * n)) * sign_eta;
      _prop_d2f[_qp] = 2.0 - n * n * (24.0 + 30.0 * n * n) * sign_eta;
      break;

    default:
      mooseError("Internal error");
  }
}