//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralPostprocessor.h"
#include "MooseVariableInterface.h"

/**
 * This postprocessor computes a volume integral of the specified variable using a range-based count.
 *
 * The integral counts areas where the absolute value of the variable falls within a specified range
 * (between low_threshold and high_threshold). Areas outside this range are excluded from the count.
 * Specializations are possible by deriving from this class and overriding computeQpIntegral().
 */
class ElementIntegralVariablePostprocessor_new1 : public ElementIntegralPostprocessor,
                                             public MooseVariableInterface<Real>
{
public:
  static InputParameters validParams();

  ElementIntegralVariablePostprocessor_new1(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  /// Holds the solution at current quadrature points
  const VariableValue & _u;
  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;
  /// The lower bound for counting (absolute value of variable must be > this)
  const Real _low_threshold;
  /// The upper bound for counting (absolute value of variable must be < this)
  const Real _high_threshold;
};