// ComputeRotatedEigenstrain.h
#pragma once

#include "ComputeEigenstrain.h"
#include "RotationTensor.h"

/**
 * ComputeRotatedEigenstrain
 *
 * Rotates a user-specified (crystal-frame) eigenstrain tensor into the global frame
 * using Bunge ZXZ Euler angles, then applies an optional scalar prefactor.
 *
 * Input params:
 *  - eigen_base (std::vector<Real>, size 6, Voigt order: xx yy zz xy xz yz; tensorial shear)
 *  - Euler_angles (RealVectorValue, default 0 0 0)  // phi1, Phi, phi2
 *  - prefactor (MaterialPropertyName or Real, default 1.0)
 */
template <bool is_ad>
class ComputeRotatedEigenstrainTempl : public ComputeEigenstrainBaseTempl<is_ad>
{
public:
  static InputParameters validParams();
  explicit ComputeRotatedEigenstrainTempl(const InputParameters & parameters);

protected:
  void computeQpEigenstrain() override;

  using ComputeEigenstrainBaseTempl<is_ad>::_qp;
  using ComputeEigenstrainBaseTempl<is_ad>::_eigenstrain;

  // Unrotated (crystal-frame) eigenstrain tensor stored locally and rotated each qp
  RankTwoTensor _eigen_base_tensor;

  // Optional scalar multiplier (can be constant or another material property)
  const GenericMaterialProperty<Real, is_ad> & _prefactor;

  // Bunge ZXZ Euler angles (phi1, Phi, phi2)
  const RealVectorValue _Euler_angles;
};

// typedefs for registration
using ComputeRotatedEigenstrain = ComputeRotatedEigenstrainTempl<false>;
using ADComputeRotatedEigenstrain = ComputeRotatedEigenstrainTempl<true>;
