//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CrystalPlasticityKalidindiUpdate_slip.h"
#include "libmesh/int_range.h"

registerMooseObject("belson324App", CrystalPlasticityKalidindiUpdate_slip);

InputParameters
CrystalPlasticityKalidindiUpdate_slip::validParams()
{
  InputParameters params = CrystalPlasticityStressUpdateBase_abs::validParams();
  params.addParam<unsigned int>(
    "slip_system_modes",
    1,
    "Number of different types of slip systems, each with its own friction strength.");
  params.addParam<std::vector<unsigned int>>(
    "number_slip_systems_per_mode",
    std::vector<unsigned int>(),
    "The number of slip systems per slip mode.");
  params.addParam<std::vector<Real>>(
    "lattice_friction_per_mode",
    std::vector<Real>(),
    "Lattice friction strength for each slip mode in MPa.");

  params.addClassDescription("Kalidindi version of homogeneous crystal plasticity.");
  params.addParam<Real>("r", 1.0, "Latent hardening coefficient");
  params.addParam<Real>("h", 541.5, "hardening constants");
  params.addParam<Real>("t_sat", 109.8, "saturated slip system strength");
  params.addParam<Real>("gss_a", 2.5, "coefficient for hardening");
  params.addParam<Real>("ao", 0.001, "slip rate coefficient");
  params.addParam<Real>("xm", 0.1, "exponent for slip rate");
  params.addParam<Real>("gss_initial", 60.8, "initial lattice friction strength of the material");

  params.addParam<MaterialPropertyName>(
      "total_twin_volume_fraction",
      "Total twin volume fraction, if twinning is considered in the simulation");

  return params;
}

CrystalPlasticityKalidindiUpdate_slip::CrystalPlasticityKalidindiUpdate_slip(
    const InputParameters & parameters)
  : CrystalPlasticityStressUpdateBase_abs(parameters),
    // Constitutive values
    _slip_system_modes(getParam<unsigned int>("slip_system_modes")),
    _number_slip_systems_per_mode(getParam<std::vector<unsigned int>>("number_slip_systems_per_mode")),
    _lattice_friction(getParam<std::vector<Real>>("lattice_friction_per_mode")),
    _r(getParam<Real>("r")),
    _h(getParam<Real>("h")),
    _tau_sat(getParam<Real>("t_sat")),
    _gss_a(getParam<Real>("gss_a")),
    _ao(getParam<Real>("ao")),
    _xm(getParam<Real>("xm")),
    _gss_initial(getParam<Real>("gss_initial")),

    // resize vectors used in the consititutive slip hardening
    _hb(_number_slip_systems, 0.0),
    _slip_resistance_increment(_number_slip_systems, 0.0),

    // resize local caching vectors used for substepping
    _previous_substep_slip_resistance(_number_slip_systems, 0.0),
    _slip_resistance_before_update(_number_slip_systems, 0.0),

    // Twinning contributions, if used
    _include_twinning_in_Lp(parameters.isParamValid("total_twin_volume_fraction")),
    _twin_volume_fraction_total(_include_twinning_in_Lp
                                    ? &getMaterialPropertyOld<Real>("total_twin_volume_fraction")
                                    : nullptr)
{
    // check that the number of slip systems is equal to the sum of the types of slip system
    if (_number_slip_systems_per_mode.size() != _slip_system_modes)
    paramError("number_slip_systems_per_mode",
               "The size the number of slip systems per mode is not equal to the number of slip "
               "system types.");                                
                                    
    if (_lattice_friction.size() != _slip_system_modes)
    paramError("lattice_friction_per_mode",
               "Please ensure that the size of lattice_friction_per_mode equals the value supplied "
               "for slip_system_modes");

  unsigned int sum = 0;
  for (const auto i : make_range(_slip_system_modes))
    sum += _number_slip_systems_per_mode[i];
  if (sum != _number_slip_systems)
    paramError("slip_system_modes",
               "The number of slip systems and the sum of the slip systems in each of the slip "
               "system modes are not equal");

}

void CrystalPlasticityKalidindiUpdate_slip::initQpStatefulProperties()
{
  CrystalPlasticityStressUpdateBase_abs::initQpStatefulProperties();
  
  // Set initial resistance from lattice friction, which is type dependent
  unsigned int slip_mode = 0;
  unsigned int counter_adjustment = 0;

  for (const auto i : make_range(_number_slip_systems))
  {
    if ((i - counter_adjustment) < _number_slip_systems_per_mode[slip_mode])
      _slip_resistance[_qp][i] = _lattice_friction[slip_mode];
    else
    {
      counter_adjustment += _number_slip_systems_per_mode[slip_mode];
      ++slip_mode;
      _slip_resistance[_qp][i] = _lattice_friction[slip_mode];
    }

    _slip_increment[_qp][i] = 0.0;
  }

}


void
CrystalPlasticityKalidindiUpdate_slip::setInitialConstitutiveVariableValues()
{
  // Would also set old dislocation densities here if included in this model
  _slip_resistance[_qp] = _slip_resistance_old[_qp];
  _previous_substep_slip_resistance = _slip_resistance_old[_qp];
}

void
CrystalPlasticityKalidindiUpdate_slip::setSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  _slip_resistance[_qp] = _previous_substep_slip_resistance;
}

bool
CrystalPlasticityKalidindiUpdate_slip::calculateSlipRate()
{
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_increment[_qp][i] =
        _ao * std::pow(std::abs(_tau[_qp][i] / _slip_resistance[_qp][i]), 1.0 / _xm);
    if (_tau[_qp][i] < 0.0)
      _slip_increment[_qp][i] *= -1.0;

    if (std::abs(_slip_increment[_qp][i]) * _substep_dt > _slip_incr_tol)
    {
      if (_print_convergence_message)
        mooseWarning("Maximum allowable slip increment exceeded ",
                     std::abs(_slip_increment[_qp][i]) * _substep_dt);

      return false;
    }
  }
  return true;
}

void
CrystalPlasticityKalidindiUpdate_slip::calculateEquivalentSlipIncrement(
    RankTwoTensor & equivalent_slip_increment)
{
  if (_include_twinning_in_Lp)
  {
    for (const auto i : make_range(_number_slip_systems))
      equivalent_slip_increment += (1.0 - (*_twin_volume_fraction_total)[_qp]) *
                                   _flow_direction[_qp][i] * _slip_increment[_qp][i] * _substep_dt;
  }
  else // if no twinning volume fraction material property supplied, use base class
    CrystalPlasticityStressUpdateBase_abs::calculateEquivalentSlipIncrement(equivalent_slip_increment);
}

void
CrystalPlasticityKalidindiUpdate_slip::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  for (const auto i : make_range(_number_slip_systems))
  {
    if (MooseUtils::absoluteFuzzyEqual(_tau[_qp][i], 0.0))
      dslip_dtau[i] = 0.0;
    else
      dslip_dtau[i] = _ao / _xm *
                      std::pow(std::abs(_tau[_qp][i] / _slip_resistance[_qp][i]), 1.0 / _xm - 1.0) /
                      _slip_resistance[_qp][i];
  }
}

bool
CrystalPlasticityKalidindiUpdate_slip::areConstitutiveStateVariablesConverged()
{
  return isConstitutiveStateVariableConverged(_slip_resistance[_qp],
                                              _slip_resistance_before_update,
                                              _previous_substep_slip_resistance,
                                              _resistance_tol);
}

void
CrystalPlasticityKalidindiUpdate_slip::updateSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  _previous_substep_slip_resistance = _slip_resistance[_qp];
}

void
CrystalPlasticityKalidindiUpdate_slip::cacheStateVariablesBeforeUpdate()
{
  _slip_resistance_before_update = _slip_resistance[_qp];
}

void CrystalPlasticityKalidindiUpdate_slip::calculateStateVariableEvolutionRateComponent()
{
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_resistance_increment[i] = 0.0;

    unsigned int slip_mode = 0;
    unsigned int counter_adjustment = 0;
    for (const auto j : make_range(_number_slip_systems))
    {
      if ((j - counter_adjustment) < _number_slip_systems_per_mode[slip_mode])
      {
        _hb[j] = _h * std::pow(std::abs(1.0 - _slip_resistance[_qp][j] / _tau_sat), _gss_a);
      }
      else
      {
        counter_adjustment += _number_slip_systems_per_mode[slip_mode];
        ++slip_mode;
        _hb[j] = _h * std::pow(std::abs(1.0 - _slip_resistance[_qp][j] / _tau_sat), _gss_a);
      }

      Real hsign = 1.0 - _slip_resistance[_qp][j] / _tau_sat;
      if (hsign < 0.0)
        _hb[j] *= -1.0;

      unsigned int iplane = i / 3;
      unsigned int jplane = j / 3;

      if (iplane == jplane) 
        _slip_resistance_increment[i] += std::abs(_slip_increment[_qp][j]) * _hb[j];
      else 
        _slip_resistance_increment[i] += std::abs(_slip_increment[_qp][j]) * _r * _hb[j];
    }
  }
}


bool
CrystalPlasticityKalidindiUpdate_slip::updateStateVariables()
{
  // Now perform the check to see if the slip system should be updated
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_resistance_increment[i] *= _substep_dt;
    if (_previous_substep_slip_resistance[i] < _zero_tol && _slip_resistance_increment[i] < 0.0)
      _slip_resistance[_qp][i] = _previous_substep_slip_resistance[i];
    else
      _slip_resistance[_qp][i] =
          _previous_substep_slip_resistance[i] + _slip_resistance_increment[i];

    if (_slip_resistance[_qp][i] < 0.0)
      return false;
  }
  return true;
}
