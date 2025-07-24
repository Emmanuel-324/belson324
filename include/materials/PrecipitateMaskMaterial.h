#pragma once

#include "Material.h"

class PrecipitateMaskMaterial : public Material
{
public:
  static InputParameters validParams();

  PrecipitateMaskMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const VariableValue & _eta1;
  const VariableValue & _eta2;
  const VariableValue & _eta3;

  MaterialProperty<Real> & _precip_mask;
};
