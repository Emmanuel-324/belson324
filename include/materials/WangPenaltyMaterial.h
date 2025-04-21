#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"

/**
 * Material that implements Wang-style phase exclusion penalty:
 * f_penalty = alpha * sum_{i<j} eta_i^2 * eta_j^2
 * This discourages spatial coexistence of multiple order parameters.
 */
class WangPenaltyMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();
  WangPenaltyMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

private:
  const unsigned int _num_eta;                                   // number of order parameters
  std::vector<const VariableValue *> _etas;                      // pointers to eta_i
  std::vector<std::string> _eta_names;                           // names of eta_i
  MaterialProperty<Real> & _f_penalty;                           // penalty free energy

  std::vector<std::vector<MaterialProperty<Real> *>> _df_deta;   // first derivatives ∂f/∂eta_i
  std::vector<std::vector<std::vector<MaterialProperty<Real> *>>> _d2f_deta2; // ∂²f/∂eta_i∂eta_j

  const Real _alpha;                                             // penalty coefficient
};
