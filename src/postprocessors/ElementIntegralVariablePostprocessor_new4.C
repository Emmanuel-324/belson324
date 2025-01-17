//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementIntegralVariablePostprocessor_new4.h"

registerMooseObject("belsonApp", ElementIntegralVariablePostprocessor_new4);

InputParameters
ElementIntegralVariablePostprocessor_new4::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");
  params.addClassDescription("Computes a volume integral of the specified variable");
  return params;
}

ElementIntegralVariablePostprocessor_new4::ElementIntegralVariablePostprocessor_new4(
    const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this,
                                 false,
                                 "variable",
                                 Moose::VarKindType::VAR_ANY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable"))
{
  addMooseVariableDependency(&mooseVariableField());
}

Real
ElementIntegralVariablePostprocessor_new4::computeQpIntegral()
{
  Real threshold = 1;  // Adjust this based on your expected \(\eta\) values for the precipitate
   return (std::abs(_u[_qp]) > threshold) ? 1.0 : 0.0; //count areas where the absolute value of eta is signifcant
}