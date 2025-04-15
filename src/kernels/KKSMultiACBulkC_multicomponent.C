//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSMultiACBulkCGeneral.h"

registerMooseObject("belson324App", KKSMultiACBulkCGeneral);


InputParameters
KKSMultiACBulkCGeneral::validParams()
{
  InputParameters params = KKSMultiACBulkBase::validParams();
  params.addClassDescription("Generalized multiphase-multicomponent kernel (part 2 of 2) for "
                             "the Bulk Allen-Cahn. Includes chemical potential terms.");
  params.addRequiredCoupledVar(
      "cj_names", "2D list of phase concentrations for each component. Each row is for one component.");
  params.addRequiredCoupledVar(
      "mu_names", "List of chemical potential variables, one for each component.");
  return params;
}

KKSMultiACBulkCGeneral::KKSMultiACBulkCGeneral(const InputParameters & parameters)
  : KKSMultiACBulkBase(parameters),
    _mu_names(coupledNames("mu_names")),
    _cj_names_2d(parameters.get<std::vector<std::vector<VariableName>>>("cj_names")),
    _n_components(_mu_names.size()),
    _prop_mu.resize(_n_components),
    _mu(_n_components),
    _cj(_n_components),
    _cj_vars(_n_components)
{
  // Sanity check
  if (_n_components != _cj_names_2d.size())
    mooseError("Number of components in mu_names must match number of cj_names rows");

  // Loop over components
  for (unsigned int a = 0; a < _n_components; ++a)
  {
    _prop_mu[a] = &getMaterialProperty<Real>(_mu_names[a]);
    _cj[a].resize(_num_j);
    _cj_vars[a].resize(_num_j);

    // Loop over phases for component a
    for (unsigned int j = 0; j < _num_j; ++j)
    {
      _cj[a][j] = &coupledValue(_cj_names_2d[a][j]);
      _cj_vars[a][j] = coupled("cj_names", a, j);
    }
  }
}

Real
KKSMultiACBulkCGeneral::computeDFDOP(PFFunctionType type)
{
  Real sum = 0.0;

  switch (type)
  {
    case Residual:
      for (unsigned int a = 0; a < _n_components; ++a)
        for (unsigned int j = 0; j < _num_j; ++j)
          sum += (*_prop_mu[a])[_qp] * (*_prop_dhjdetai[j])[_qp] * (*_cj[a][j])[_qp];

      return -sum;

    case Jacobian:
      // Only contribute if eta_i is the nonlinear variable
      if (_etai_var != _var.number())
        return 0.0;

      for (unsigned int a = 0; a < _n_components; ++a)
        for (unsigned int j = 0; j < _num_j; ++j)
          sum += (*_prop_mu[a])[_qp] * (*_prop_d2hjdetai2[j])[_qp] * (*_cj[a][j])[_qp];

      return -_phi[_j][_qp] * sum;
  }

  mooseError("Invalid type passed in");
}

Real
KKSMultiACBulkCGeneral::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real res = ACBulk<Real>::computeQpOffDiagJacobian(jvar);
  Real sum = 0.0;

  // Loop over all components and phases
  for (unsigned int a = 0; a < _n_components; ++a)
  {
    for (unsigned int j = 0; j < _num_j; ++j)
    {
      if (jvar == _cj_vars[a][j])
      {
        // Derivative w.r.t. cj[a][j]
        res -= _L[_qp] * (*_prop_mu[a])[_qp] * (*_prop_dhjdetai[j])[_qp] *
               _phi[_j][_qp] * _test[_i][_qp];
        return res;
      }
    }
  }

  // For all other vars, treat as general arguments (e.g., temperature, other phase fields)
  const unsigned int cvar = mapJvarToCvar(jvar);

  for (unsigned int a = 0; a < _n_components; ++a)
  {
    for (unsigned int j = 0; j < _num_j; ++j)
    {
      sum += (*_prop_mu[a])[_qp] * (*_prop_d2hjdetaidarg[j][cvar])[_qp] * (*_cj[a][j])[_qp];
    }
  }

  res -= _L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];

  return res;
}

