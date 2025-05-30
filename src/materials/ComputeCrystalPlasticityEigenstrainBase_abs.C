//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeCrystalPlasticityEigenstrainBase_abs.h"
#include "RankTwoTensor.h"

InputParameters
ComputeCrystalPlasticityEigenstrainBase_abs::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();

  // The return stress increment classes are intended to be iterative materials, so must set compute
  // = false for all inheriting classes
  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");

  params.addRequiredParam<std::string>(
      "deformation_gradient_name",
      "Material property name for the deformation gradient tensor computed "
      "by this model.");
  return params;
}

ComputeCrystalPlasticityEigenstrainBase_abs::ComputeCrystalPlasticityEigenstrainBase_abs(
    const InputParameters & parameters)
  : ComputeEigenstrainBase(parameters),
    _deformation_gradient_name(_base_name + getParam<std::string>("deformation_gradient_name")),
    _deformation_gradient(declareProperty<RankTwoTensor>(_deformation_gradient_name)),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>(_deformation_gradient_name)),
    _crysrot(getMaterialProperty<RankTwoTensor>(
      _base_name +   "crysrot")) // defined in the elasticity tensor classes for crystal plasticity
{
}

void
ComputeCrystalPlasticityEigenstrainBase_abs::initQpStatefulProperties()
{
  ComputeEigenstrainBase::initQpStatefulProperties();

  // Initialize deformation gradient to be identity
  _deformation_gradient[_qp].zero();
  _deformation_gradient[_qp].addIa(1.0);
}

void
ComputeCrystalPlasticityEigenstrainBase_abs::computeQpEigenstrain()
{
  // updates the deformation gradient from the child class
  computeQpDeformationGradient();
  // use the updated deformation gradient for the eigenstrain calculation
  _eigenstrain[_qp] = 0.5 * (_deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] -
                             RankTwoTensor::Identity());
}

void
ComputeCrystalPlasticityEigenstrainBase_abs::setQp(const unsigned int & qp)
{
  _qp = qp;
}

void
ComputeCrystalPlasticityEigenstrainBase_abs::setSubstepDt(const Real & substep_dt)
{
  _substep_dt = substep_dt;
}

const RankTwoTensor
ComputeCrystalPlasticityEigenstrainBase_abs::getDeformationGradient() const
{
  return _deformation_gradient[_qp];
}

const RankTwoTensor
ComputeCrystalPlasticityEigenstrainBase_abs::getDeformationGradientInverse() const
{
  return _deformation_gradient[_qp].inverse();
}
