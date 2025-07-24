#include "PrecipitateMaskMaterial.h"

registerMooseObject("belson324App", PrecipitateMaskMaterial);

InputParameters
PrecipitateMaskMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("eta1", "Order parameter for precipitate variant 1");
  params.addRequiredCoupledVar("eta2", "Order parameter for precipitate variant 2");
  params.addRequiredCoupledVar("eta3", "Order parameter for precipitate variant 3");
  params.addClassDescription("Computes total precipitate area mask as eta1^2 + eta2^2 + eta3^2");
  return params;
}

PrecipitateMaskMaterial::PrecipitateMaskMaterial(const InputParameters & parameters)
  : Material(parameters),
    _eta1(coupledValue("eta1")),
    _eta2(coupledValue("eta2")),
    _eta3(coupledValue("eta3")),
    _precip_mask(declareProperty<Real>("precip_mask"))
{
}


void
PrecipitateMaskMaterial::computeQpProperties()
{
  _precip_mask[_qp] = _eta1[_qp] * _eta1[_qp]
                    + _eta2[_qp] * _eta2[_qp]
                    + _eta3[_qp] * _eta3[_qp];
}
