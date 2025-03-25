//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MatReaction_abscouple.h"

registerMooseObject("belson324App", MatReaction_abscouple);


InputParameters
MatReaction_abscouple::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addCoupledVar("v",
                       "Set this to make v a coupled variable, otherwise it will use the "
                       "kernel's nonlinear variable for v");
  params.addClassDescription("Kernel to add -L*v, where L=reaction rate, v=variable");
  params.addParam<MaterialPropertyName>("reaction_rate", "L", "The reaction rate used with the kernel");
  params.addCoupledVar("args", "Vector of nonlinear variable arguments this object depends on");
  return params;
}

MatReaction_abscouple::MatReaction_abscouple(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    _is_coupled(isCoupled("v")),
    _v_name(_is_coupled ? getVar("v", 0)->name() : _var.name()),
    _v(_is_coupled ? coupledValue("v") : _u),
    _v_var(_is_coupled ? coupled("v") : _var.number()),
    _L(getMaterialProperty<Real>("reaction_rate")),
    _eta_name(_var.name()),
    _dLdop(getMaterialPropertyDerivative<Real>("reaction_rate", _eta_name)),
    _dLdv(getMaterialPropertyDerivative<Real>("reaction_rate", _v_name)),
    _nvar(_coupled_moose_vars.size()),
    _dLdarg(_nvar)
{
  // Get reaction rate derivatives
  for (unsigned int i = 0; i < _nvar; ++i)
  {
    MooseVariableFEBase * ivar = _coupled_moose_vars[i];
    _dLdarg[i] = &getMaterialPropertyDerivative<Real>("reaction_rate", ivar->name());
  }
}

void
MatReaction_abscouple::initialSetup()
{
  validateNonlinearCoupling<Real>("reaction_rate");
}

Real
MatReaction_abscouple::computeQpResidual()
{
  int sign_v = 0;
  if (_v[_qp] > 0) sign_v = 1;
  else if (_v[_qp] < 0) sign_v = -1;

  return -_L[_qp] * _test[_i][_qp] * _v[_qp] * sign_v;
}

Real
MatReaction_abscouple::computeQpJacobian()
{
  int sign_v = 0;
  if (_v[_qp] > 0) sign_v = 1;
  else if (_v[_qp] < 0) sign_v = -1;

  if (_is_coupled)
    return -_dLdop[_qp] * _v[_qp] * _phi[_j][_qp] * _test[_i][_qp] * sign_v;

  return -(_L[_qp] + _dLdop[_qp] * _v[_qp]) * _phi[_j][_qp] * _test[_i][_qp] * sign_v;
}

Real
MatReaction_abscouple::computeQpOffDiagJacobian(unsigned int jvar)
{
  // first handle the case where jvar is a coupled variable v being added to residual
  // the first term in the sum just multiplies by L which is always needed
  // the second term accounts for cases where L depends on v
  int sign_v = 0;
  if (_v[_qp] > 0) sign_v = 1;
  else if (_v[_qp] < 0) sign_v = -1;

  if ((jvar == _v_var) && _is_coupled)
    return -(_L[_qp] + _dLdv[_qp] * _v[_qp]) * _test[_i][_qp] * _phi[_j][_qp] * sign_v;

  //  for all other vars get the coupled variable jvar is referring to
  const unsigned int cvar = mapJvarToCvar(jvar);

  return -(*_dLdarg[cvar])[_qp] * _v[_qp] * _test[_i][_qp] * _phi[_j][_qp] * sign_v;
}
