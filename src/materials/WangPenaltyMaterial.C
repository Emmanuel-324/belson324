// WangPenaltyMaterial.C
#include "WangPenaltyMaterial.h"

registerMooseObject("belson324App", WangPenaltyMaterial);

InputParameters
WangPenaltyMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Wang-style penalty term to discourage phase coexistence: alpha * sum_{i<j} eta_i^2 * eta_j^2");
  params.addRequiredCoupledVar("etas", "List of order parameters (eta_i)");
  params.addParam<Real>("alpha", 1000.0, "Penalty prefactor alpha");
  return params;
}

WangPenaltyMaterial::WangPenaltyMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _num_eta(coupledComponents("etas")),
    _etas(_num_eta),
    _eta_names(coupledNames("etas")),
    _f_penalty(declareProperty<Real>("f_penalty")),
    _df_deta(_num_eta),
    _d2f_deta2(_num_eta),
    _alpha(getParam<Real>("alpha"))
{
  for (unsigned int i = 0; i < _num_eta; ++i)
  {
    _etas[i] = &coupledValue("etas", i);
    _df_deta[i].resize(_num_eta);
    _d2f_deta2[i].resize(_num_eta);
    for (unsigned int j = 0; j < _num_eta; ++j)
    {
      _df_deta[i][j] = &declarePropertyDerivative<Real>("f_penalty", _eta_names[i]);
      _d2f_deta2[i][j].resize(_num_eta);
      for (unsigned int k = 0; k < _num_eta; ++k)
        _d2f_deta2[i][j][k] = &declarePropertyDerivative<Real>("f_penalty", _eta_names[i], _eta_names[k]);
    }
  }
}

void WangPenaltyMaterial::computeQpProperties()
{
  Real penalty = 0.0;

  // Compute penalty energy
  for (unsigned int i = 0; i < _num_eta; ++i)
    for (unsigned int j = i + 1; j < _num_eta; ++j)
      penalty += (*_etas[i])[_qp] * (*_etas[i])[_qp] * (*_etas[j])[_qp] * (*_etas[j])[_qp];

  _f_penalty[_qp] = _alpha * penalty;

  // Compute derivatives
  for (unsigned int i = 0; i < _num_eta; ++i)
  {
    Real value_i = (*_etas[i])[_qp];
    for (unsigned int j = 0; j < _num_eta; ++j)
    {
      Real value_j = (*_etas[j])[_qp];
      if (i != j)
        (*_df_deta[i][j])[_qp] = 2.0 * _alpha * value_i * value_j * value_j;
      else
      {
        Real sum = 0.0;
        for (unsigned int k = 0; k < _num_eta; ++k)
          if (k != i)
            sum += value_k * value_k;
        (*_df_deta[i][j])[_qp] = 2.0 * _alpha * value_i * sum;
      }

      // Second derivatives (simplified symmetrical version)
      for (unsigned int k = 0; k < _num_eta; ++k)
        (*_d2f_deta2[i][j][k])[_qp] = 0.0; // optionally populate if needed
    }
  }
}
