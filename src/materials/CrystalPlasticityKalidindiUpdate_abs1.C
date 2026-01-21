//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html



#include "CrystalPlasticityKalidindiUpdate_abs1.h"
#include "libmesh/int_range.h"

registerMooseObject("belson324App", CrystalPlasticityKalidindiUpdate_abs1);

InputParameters
CrystalPlasticityKalidindiUpdate_abs1::validParams()
{
  InputParameters params = CrystalPlasticityStressUpdateBase_abs1::validParams();
  params.addClassDescription("Kalidindi version of homogeneous crystal plasticity.");
  params.addParam<Real>("r", 1.0, "Latent hardening coefficient");
  params.addParam<Real>("h", 541.5, "hardening constants");
  params.addParam<Real>("t_sat", 109.8, "saturated slip system strength");
  params.addParam<Real>("gss_a", 2.5, "coefficient for hardening");
  params.addParam<Real>("ao", 0.001, "slip rate coefficient");
  params.addParam<Real>("xm", 0.1, "exponent for slip rate");
  params.addParam<Real>("gss_initial", 60.8, "initial lattice friction strength of the material");
  params.addParam<Real>("c_back", 1000.0, "Backstress hardening modulus");
  params.addParam<Real>("d_back", 10.0, "Backstress dynamic recovery");
  params.addParam<Real>("backstress_tol", 1e-6, "Backstress convergence tolerance");



  params.addParam<MaterialPropertyName>(
      "total_twin_volume_fraction",
      "Total twin volume fraction, if twinning is considered in the simulation");

  return params;
}

CrystalPlasticityKalidindiUpdate_abs1::CrystalPlasticityKalidindiUpdate_abs1(
    const InputParameters & parameters)
  : CrystalPlasticityStressUpdateBase_abs1(parameters),
    // Constitutive values
    _r(getParam<Real>("r")),
    _h(getParam<Real>("h")),
    _tau_sat(getParam<Real>("t_sat")),
    _gss_a(getParam<Real>("gss_a")),
    _ao(getParam<Real>("ao")),
    _xm(getParam<Real>("xm")),
    _gss_initial(getParam<Real>("gss_initial")),
    _c_back(getParam<Real>("c_back")),
    _d_back(getParam<Real>("d_back")),
    _backstress_tol(getParam<Real>("backstress_tol")),

    // resize vectors used in the consititutive slip hardening
    _hb(_number_slip_systems, 0.0),
    _slip_resistance_increment(_number_slip_systems, 0.0),

    // resize local caching vectors used for substepping
    _previous_substep_slip_resistance(_number_slip_systems, 0.0),
    _slip_resistance_before_update(_number_slip_systems, 0.0),
    _previous_substep_backstress(_number_slip_systems, 0.0),
    _backstress_before_update(_number_slip_systems, 0.0),
    _backstress_increment(_number_slip_systems, 0.0),


    // Twinning contributions, if used
    _include_twinning_in_Lp(parameters.isParamValid("total_twin_volume_fraction")),
    _twin_volume_fraction_total(_include_twinning_in_Lp
                                    ? &getMaterialPropertyOld<Real>("total_twin_volume_fraction")
                                    : nullptr)
{
}

void
CrystalPlasticityKalidindiUpdate_abs1::initQpStatefulProperties()
{
  CrystalPlasticityStressUpdateBase_abs1::initQpStatefulProperties();
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_resistance[_qp][i] = _gss_initial;
    _slip_increment[_qp][i] = 0.0;
  }
}

void
CrystalPlasticityKalidindiUpdate_abs1::setInitialConstitutiveVariableValues()
{
  // Would also set old dislocation densities here if included in this model
  _slip_resistance[_qp] = _slip_resistance_old[_qp];
  _previous_substep_slip_resistance = _slip_resistance_old[_qp];

  // NEW: mirror slip_resistance behavior for backstress
  _backstress[_qp] = _backstress_old[_qp];
  _previous_substep_backstress = _backstress_old[_qp];
}

void
CrystalPlasticityKalidindiUpdate_abs1::setSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  _slip_resistance[_qp] = _previous_substep_slip_resistance;
  // NEW: mirror slip_resistance behavior for backstress
  _backstress[_qp] = _previous_substep_backstress;
}

bool
CrystalPlasticityKalidindiUpdate_abs1::calculateSlipRate()
{
  for (const auto i : make_range(_number_slip_systems))
  {
    const Real tau_eff = _tau[_qp][i] - _backstress[_qp][i];
    _slip_increment[_qp][i] =
        _ao * std::pow(std::abs(tau_eff / _slip_resistance[_qp][i]), 1.0 / _xm);
    if (tau_eff < 0.0)
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
CrystalPlasticityKalidindiUpdate_abs1::calculateEquivalentSlipIncrement(
    RankTwoTensor & equivalent_slip_increment)
{
  if (_include_twinning_in_Lp)
  {
    for (const auto i : make_range(_number_slip_systems))
      equivalent_slip_increment += (1.0 - (*_twin_volume_fraction_total)[_qp]) *
                                   _flow_direction[_qp][i] * _slip_increment[_qp][i] * _substep_dt;
  }
  else // if no twinning volume fraction material property supplied, use base class
    CrystalPlasticityStressUpdateBase_abs1::calculateEquivalentSlipIncrement(equivalent_slip_increment);
}

void
CrystalPlasticityKalidindiUpdate_abs1::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  for (const auto i : make_range(_number_slip_systems))
  {
    const Real tau_eff = _tau[_qp][i] - _backstress[_qp][i];

    if (MooseUtils::absoluteFuzzyEqual(tau_eff, 0.0))
      dslip_dtau[i] = 0.0;
    else
      dslip_dtau[i] = _ao / _xm *
          std::pow(std::abs(tau_eff / _slip_resistance[_qp][i]), 1.0 / _xm - 1.0) /
          _slip_resistance[_qp][i];
  }
}

bool
CrystalPlasticityKalidindiUpdate_abs1::areConstitutiveStateVariablesConverged()
{
  const bool slip_ok =
      isConstitutiveStateVariableConverged(_slip_resistance[_qp],
                                           _slip_resistance_before_update,
                                           _previous_substep_slip_resistance,
                                           _resistance_tol);

  const bool back_ok =
      isConstitutiveStateVariableConverged(_backstress[_qp],
                                           _backstress_before_update,
                                           _previous_substep_backstress,
                                           _backstress_tol);

  return slip_ok && back_ok;
}


void
CrystalPlasticityKalidindiUpdate_abs1::updateSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  _previous_substep_slip_resistance = _slip_resistance[_qp];
  // NEW: mirror slip_resistance behavior for backstress
  _previous_substep_backstress = _backstress[_qp];
}

void
CrystalPlasticityKalidindiUpdate_abs1::cacheStateVariablesBeforeUpdate()
{
  _slip_resistance_before_update = _slip_resistance[_qp];
  // NEW: mirror slip_resistance behavior for backstress
  _backstress_before_update = _backstress[_qp];
}

void
CrystalPlasticityKalidindiUpdate_abs1::calculateStateVariableEvolutionRateComponent()
{
  for (const auto i : make_range(_number_slip_systems))
  {  
     const Real gdot = _slip_increment[_qp][i]; // this is your slip rate
    _backstress_increment[i] = _c_back * gdot - _d_back * std::abs(gdot) * _backstress[_qp][i];
    // Clear out increment from the previous iteration
    _slip_resistance_increment[i] = 0.0;

    _hb[i] = _h * std::pow(std::abs(1.0 - _slip_resistance[_qp][i] / _tau_sat), _gss_a);
    const Real hsign = 1.0 - _slip_resistance[_qp][i] / _tau_sat;
    if (hsign < 0.0)
      _hb[i] *= -1.0;
  }

  for (const auto i : make_range(_number_slip_systems))
  {
    for (const auto j : make_range(_number_slip_systems))
    {
      unsigned int iplane, jplane;
      iplane = i / 3;
      jplane = j / 3;

      if (iplane == jplane) // self vs. latent hardening
        _slip_resistance_increment[i] +=
            std::abs(_slip_increment[_qp][j]) * _hb[j]; // q_{ab} = 1.0 for self hardening
      else
        _slip_resistance_increment[i] +=
            std::abs(_slip_increment[_qp][j]) * _r * _hb[j]; // latent hardenign
    }
  }
}

bool
CrystalPlasticityKalidindiUpdate_abs1::updateStateVariables()
{
  // Now perform the check to see if the slip system should be updated
  for (const auto i : make_range(_number_slip_systems))
  {
    _backstress[_qp][i] =
    _previous_substep_backstress[i] + _backstress_increment[i] * _substep_dt;
    _slip_resistance_increment[i] *= _substep_dt;
    if (_previous_substep_slip_resistance[i] < _zero_tol && _slip_resistance_increment[i] < 0.0)
      _slip_resistance[_qp][i] = _previous_substep_slip_resistance[i];
    else
      _slip_resistance[_qp][i] =
          _previous_substep_slip_resistance[i] + _slip_resistance_increment[i];

    if (_slip_resistance[_qp][i] < 0.0)
      return false;
    if (!std::isfinite(_backstress[_qp][i]))
      return false;
  }
  return true;
}
