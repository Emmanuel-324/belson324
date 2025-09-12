#include "ComputeRotatedEigenstrain.h"
#include "RotationTensor.h"

registerMooseObject("belson324App", ComputeRotatedEigenstrain);
registerMooseObject("belson324App", ADComputeRotatedEigenstrain);

template <bool is_ad>
InputParameters
ComputeRotatedEigenstrainTempl<is_ad>::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();
  params.addClassDescription("Computes an eigenstrain rotated by Euler angles");
  params.addRequiredParam<std::vector<Real>>(
      "eigen_base",
      "Unrotated (crystal-frame) eigenstrain, Voigt order: xx yy zz xy xz yz (tensorial shear)");
  params.addParam<RealVectorValue>(
      "Euler_angles", RealVectorValue(0, 0, 0),
      "Bunge ZXZ Euler angles (phi1, Phi, phi2) to rotate eigenstrain into global frame");
  params.addParam<MaterialPropertyName>(
      "prefactor", 1.0, "Material property or constant scaling the rotated eigenstrain");
  return params;
}

template <bool is_ad>
ComputeRotatedEigenstrainTempl<is_ad>::ComputeRotatedEigenstrainTempl(
    const InputParameters & parameters)
  : ComputeEigenstrainBaseTempl<is_ad>(parameters),
    _prefactor(this->template getGenericMaterialProperty<Real, is_ad>("prefactor")),
    _Euler_angles(this->template getParam<RealVectorValue>("Euler_angles"))
{
  _eigen_base_tensor.fillFromInputVector(
      this->template getParam<std::vector<Real>>("eigen_base"));
}

template <bool is_ad>
void
ComputeRotatedEigenstrainTempl<is_ad>::computeQpEigenstrain()
{
  RotationTensor R(_Euler_angles);              // Bunge ZXZ rotation
  RankTwoTensor rotated = R * _eigen_base_tensor * R.transpose();
  this->_eigenstrain[_qp] = rotated * _prefactor[_qp];
}

template class ComputeRotatedEigenstrainTempl<false>;
template class ComputeRotatedEigenstrainTempl<true>;
