//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSMultiPhaseConcentrationTernary.h"

registerMooseObject("belson324App", KKSMultiPhaseConcentrationTernary);

InputParameters
KKSMultiPhaseConcentrationTernary::validParams()
{
  InputParameters params = KernelValue::validParams();
  params.addClassDescription(
      "KKS multi-phase model kernel for ternary alloys enforcing concentration interpolation for two solutes.");

  params.addRequiredCoupledVar("c_Al", "Physical concentration of Al");
  params.addRequiredCoupledVar("c_Nb", "Physical concentration of Nb");

  params.addRequiredCoupledVar(
      "cj_Al", "Array of Al phase concentrations cj_Al. Place in same order as hj_Al_names!");
  params.addRequiredCoupledVar(
      "cj_Nb", "Array of Nb phase concentrations cj_Nb. Place in same order as hj_Nb_names!");

  params.addCoupledVar("etas", "Order parameters for all phases");

  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "hj_Al_names", "Switching Function Materials for Al concentration.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "hj_Nb_names", "Switching Function Materials for Nb concentration.");

  return params;
}

KKSMultiPhaseConcentrationTernary::KKSMultiPhaseConcentrationTernary(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<KernelValue>>(parameters),
    _num_j(coupledComponents("cj_Al")),
    _cj_Al(coupledValues("cj_Al")),
    _cj_Nb(coupledValues("cj_Nb")),
    _cj_Al_map(getParameterJvarMap("cj_Al")),  
    _cj_Nb_map(getParameterJvarMap("cj_Nb")),  
    _eta_map(getParameterJvarMap("etas")),    
    _c_Al(coupledValue("c_Al")),
    _c_Nb(coupledValue("c_Nb")),
    _hj_Al_names(getParam<std::vector<MaterialPropertyName>>("hj_Al_names")),
    _hj_Nb_names(getParam<std::vector<MaterialPropertyName>>("hj_Nb_names")),
    _prop_hj_Al(_hj_Al_names.size()),
    _prop_hj_Nb(_hj_Nb_names.size()),
    _k(-1) 
{
  // Ensure consistency in number of phases
  if (_num_j != _hj_Al_names.size() || _num_j != _hj_Nb_names.size())
    mooseError("Mismatch in number of phases for concentration interpolation");

  // Load material properties for switching functions
  for (unsigned int m = 0; m < _num_j; ++m)
  {
    _prop_hj_Al[m] = &getMaterialPropertyByName<Real>(_hj_Al_names[m]);
    _prop_hj_Nb[m] = &getMaterialPropertyByName<Real>(_hj_Nb_names[m]);
  }
  //  Find which phase concentration this kernel is solving for
  for (unsigned int m = 0; m < _num_j; ++m)
  {
    if (coupled("cj_Al", m) == _var.number())
    {
      _k = m;
      break;  // Stop once we find the right variable
    }
  }

  // Error if no valid nonlinear variable is found
  if (_k < 0)
    mooseError("Nonlinear variable must be one of the phase concentrations (cj_Al or cj_Nb).");
}

Real
KKSMultiPhaseConcentrationTernary::precomputeQpResidual()
{
  // Residuals for Al and Nb concentrations
  Real sum_ch_Al = 0.0, sum_ch_Nb = 0.0;

  for (unsigned int m = 0; m < _num_j; ++m)
  {
    sum_ch_Al += (*_cj_Al[m])[_qp] * (*_prop_hj_Al[m])[_qp];
    sum_ch_Nb += (*_cj_Nb[m])[_qp] * (*_prop_hj_Nb[m])[_qp];
  }

  return sum_ch_Al - _c_Al[_qp] + sum_ch_Nb - _c_Nb[_qp];
}

Real
KKSMultiPhaseConcentrationTernary::precomputeQpJacobian()
{
  return (*_prop_hj_Al[_k])[_qp] * _phi[_j][_qp] + (*_prop_hj_Nb[_k])[_qp] * _phi[_j][_qp];
}

Real
KKSMultiPhaseConcentrationTernary::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == coupled("c_Al"))
    return -_test[_i][_qp] * _phi[_j][_qp];

  else if (jvar == coupled("c_Nb"))
    return -_test[_i][_qp] * _phi[_j][_qp];

  auto cjvar_Al = mapJvarToCvar(jvar, _cj_Al_map);
  if (cjvar_Al >= 0)
    return _test[_i][_qp] * (*_prop_hj_Al[cjvar_Al])[_qp] * _phi[_j][_qp];

  auto cjvar_Nb = mapJvarToCvar(jvar, _cj_Nb_map);
  if (cjvar_Nb >= 0)
    return _test[_i][_qp] * (*_prop_hj_Nb[cjvar_Nb])[_qp] * _phi[_j][_qp];

  return 0.0;
}

