//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "DelimitedFileReader.h"

class CrystalPlasticityStressUpdateBase_new : public CrystalPlasticityStressUpdateBase
{
public:
  static InputParameters validParams();

  CrystalPlasticityStressUpdateBase_new(const InputParameters & parameters);

  /// Declare slip system plane normals and directions as static members
  static const std::vector<std::array<double, 3>> slip_planes;
  static const std::vector<std::array<double, 3>> slip_directions;

  /// Sets the value of the global variable _qp for inheriting classes
  void setQp(const unsigned int & qp);

  /// Sets the value of the _substep_dt for inheriting classes
  void setSubstepDt(const Real & substep_dt);

  ///@{ Retained as empty methods to avoid a warning from Material.C in framework. These methods are unused in all inheriting classes and should not be overwritten.
  virtual void resetQpProperties() final {}
  virtual void resetProperties() final {}
  ///@}

  /**
   * Initializes the stateful properties such as PK2 stress, resolved shear
   * stress, plastic deformation gradient, slip system resistances, etc.
   * This class is often overwritten by inheriting classes.
   */
  virtual void initQpStatefulProperties() override;
  virtual void setMaterialVectorSize();

  /**
   * A helper method to read in plane normal and direction vectors from a file
   * and to normalize the vectors. This method is abstracted to allow for reuse
   * in inheriting classes with multiple plane normal and direction vector pairs.
   */
  virtual void getSlipSystems();

  /**
   * A helper method to transform the Miller-Bravais 4-index notation for HCP
   * crystals into a 3-index Cartesian representation, using the convention
   * a$_1$ = x of axis alignment in the basal plane.
   */
  void transformHexagonalMillerBravaisSlipSystems(const MooseUtils::DelimitedFileReader & reader);

  /**
   * Sorts slip systems specifically for DO22 (Tetragonal) crystal structures,
   * categorizing slip systems into primary and secondary sets.
   */
  void sortSlipSystemsForDO22();

  /**
   * Computes the Schmid tensor (m x n) for the original (reference) crystal
   * lattice orientation for each glide slip system.
   */
  void calculateFlowDirection(const RankTwoTensor & crysrot);

  /**
   * Computes the shear stress for each slip system.
   */
  void calculateShearStress(const RankTwoTensor & pk2,
                            const RankTwoTensor & inverse_eigenstrain_deformation_grad,
                            const unsigned int & num_eigenstrains);

  /**
   * Calculates the total value of $\frac{d\mathbf{F}^P^{-1}}{d\mathbf{PK2}}$ and
   * is intended to be an overwritten helper method for inheriting classes with
   * multiple constitutive dislocation slip mechanisms, e.g. glide and twinning.
   */
  virtual void calculateTotalPlasticDeformationGradientDerivative(
      RankFourTensor & dfpinvdpk2,
      const RankTwoTensor & inverse_plastic_deformation_grad_old,
      const RankTwoTensor & inverse_eigenstrain_deformation_grad_old,
      const unsigned int & num_eigenstrains);

  /**
   * A helper method to rotate the slip direction and plane normal into
   * the local crystal lattice orientation as defined by the crystal rotation
   * tensor from the Elasticity tensor class.
   */
  void calculateSchmidTensor(const unsigned int & number_dislocation_systems,
                             const std::vector<RealVectorValue> & plane_normal_vector,
                             const std::vector<RealVectorValue> & direction_vector,
                             std::vector<RankTwoTensor> & schmid_tensor,
                             const RankTwoTensor & crysrot);

  /**
   * A helper method to sort the slip systems of a crystal into cross slip families based
   * on common slip directions. This method determines if slip directions are parallel
   * and stores the index of the slip systems in a vector of vectors, where the outer vector
   * separates the individual slip system families.
   */
  void sortCrossSlipFamilies();

  /**
   * Identifies to which cross-slip family a particular slip system index belongs
   * after the slip systems have been sorted. Returns the integer value of the
   * identified cross-slip system family.
   */
  unsigned int identifyCrossSlipFamily(const unsigned int index);

  /**
   * Sets the constitutive internal state variables' current value and the previous substep
   * value to the old property value for the start of the stress convergence while loop.
   */
  virtual void setInitialConstitutiveVariableValues() {}

  /**
   * Sets the current constitutive internal state variable value to that of the
   * previous substep at the beginning of the next substep increment.
   */
  virtual void setSubstepConstitutiveVariableValues() {}

  /**
   * Stores the current value of the constitutive internal state variables into
   * a separate material property in case substepping is required.
   */
  virtual void updateSubstepConstitutiveVariableValues() {}

  /**
   * This virtual method is called to calculate the slip system slip
   * increment based on the constitutive model defined in the child class.
   * This method must be overwritten in the child class.
   */
  virtual bool calculateSlipRate() = 0;

  /**
   * Computes the equivalent slip increment across all slip systems.
   */
  virtual void calculateEquivalentSlipIncrement(RankTwoTensor & /*equivalent_slip_increment*/);

  /**
   * Finds the derivative of the slip increment with respect to the applied shear
   * stress on the slip system based on the constitutive model defined in the child class.
   */
  virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & /*dslip_dtau*/) = 0;

  /**
   * Finalizes the values of the state variables and slip system resistance
   * for the current timestep after convergence has been reached.
   */
  virtual void cacheStateVariablesBeforeUpdate() {}
  virtual void calculateStateVariableEvolutionRateComponent() {}
  virtual bool updateStateVariables() = 0;
  virtual void calculateSlipResistance() {}

  /**
   * Determines if all the state variables have converged.
   */
  virtual bool areConstitutiveStateVariablesConverged() { return true; }

  /**
   * Checks if a typical state variable, e.g., defect density, has converged
   * by comparing the change in the values over the iteration period.
   */
  virtual bool isConstitutiveStateVariableConverged(const std::vector<Real> & current_var,
                                                    const std::vector<Real> & var_before_update,
                                                    const std::vector<Real> & previous_substep_var,
                                                    const Real & tolerance);

protected:
  /// Base name prepended to all material property names to allow for multi-material systems
  const std::string _base_name;

  /// **Updated to include TETRAGONAL (DO22)**
  const enum class CrystalLatticeType { BCC, FCC, HCP, TETRAGONAL } _crystal_lattice_type;
  
  std::vector<unsigned int> _sorted_slip_systems;  // Stores sorted slip system indices

  const std::vector<Real> _unit_cell_dimension;
  const unsigned int _number_slip_systems;
  std::string _slip_sys_file_name;

  const Real _number_cross_slip_directions;
  const Real _number_cross_slip_planes;

  Real _rel_state_var_tol;
  Real _slip_incr_tol;
  Real _resistance_tol;
  Real _zero_tol;

  MaterialProperty<std::vector<Real>> & _slip_resistance;
  const MaterialProperty<std::vector<Real>> & _slip_resistance_old;
  MaterialProperty<std::vector<Real>> & _slip_increment;

  std::vector<RealVectorValue> _slip_direction;
  std::vector<RealVectorValue> _slip_plane_normal;
  MaterialProperty<std::vector<RankTwoTensor>> & _flow_direction;

  MaterialProperty<std::vector<Real>> & _tau;
  const bool _print_convergence_message;
  Real _substep_dt;

  std::vector<std::vector<unsigned int>> _cross_slip_familes;
  bool _calculate_cross_slip;
};
