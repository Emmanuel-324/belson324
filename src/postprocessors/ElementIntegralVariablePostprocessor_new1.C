//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementIntegralVariablePostprocessor_new1.h"

registerMooseObject("belson324App", ElementIntegralVariablePostprocessor_new1);

InputParameters
ElementIntegralVariablePostprocessor_new1::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");
  params.addParam<Real>("low_threshold", 0.5, "The lower bound for counting (absolute value of variable must be > this)");
  params.addParam<Real>("high_threshold", 1.0, "The upper bound for counting (absolute value of variable must be < this)");
  params.addClassDescription("Computes a volume integral of the specified variable using a range-based count");
  return params;
}

ElementIntegralVariablePostprocessor_new1::ElementIntegralVariablePostprocessor_new1(
    const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this,
                                 false,
                                 "variable",
                                 Moose::VarKindType::VAR_ANY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable")),
    _low_threshold(getParam<Real>("low_threshold")),
    _high_threshold(getParam<Real>("high_threshold"))
{
  addMooseVariableDependency(&mooseVariableField());
}

Real
ElementIntegralVariablePostprocessor_new1::computeQpIntegral()
{
  Real abs_u = std::abs(_u[_qp]);
  return (abs_u > _low_threshold && abs_u < _high_threshold) ? 1.0 : 0.0; // Count only within the range
}