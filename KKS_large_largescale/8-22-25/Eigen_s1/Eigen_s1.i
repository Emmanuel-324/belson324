# This test is for the multicomponent In718 alloy
[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [phasem_gr0]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 175
    ny = 175
#   nz = 2
    xmin = 0
    xmax = 350
    ymin = 0
    ymax = 350
    zmin = 0
    zmax = 0
    elem_type = QUAD4
  []
  [phasem_gr0_id]
    type = SubdomainIDGenerator
    input = phasem_gr0
    subdomain_id = 0
  []
  [phasem_gr1]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 175
    ny = 175
#   nz = 2
    xmin = 350
    xmax = 700
    ymin = 0
    ymax = 350
    zmin = 0
    zmax = 0
    elem_type = QUAD4
  []
  [phasem_gr1_id]
    type = SubdomainIDGenerator
    input = phasem_gr1
    subdomain_id = 1
  []
  [sticher]
    type = StitchedMeshGenerator
    inputs = 'phasem_gr0_id phasem_gr1_id'
    stitch_boundaries_pairs = 'right left'
    prevent_boundary_ids_overlap = false
  []
[]

[BCs]
  [./u_right_pull]
    type = Pressure
    displacements = 'disp_x disp_y'
    variable = disp_x
    boundary = right
    factor = 0
  [../]
  [./all_c1]
    type =  NeumannBC
    variable = 'c1'
    boundary = 'left right top bottom'
    value = 0
  [../]	
  [./all_c2]
    type =  NeumannBC
    variable = 'c2'
    boundary = 'left right top bottom'
    value = 0
  [../]	
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
[]


[AuxVariables]
  [./gb_scale_aux]
    family = MONOMIAL
    order  = FIRST
  [../]
  [./Energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx_gr0]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./stress_xx_gr1]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./e_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]
  [temperature_gr0]
    order = FIRST
    family = LAGRANGE
    block = 0 
  []
  [temperature_gr1]
    order = FIRST
    family = LAGRANGE
    block = 1
  []
[]

[Bounds]
  [./eta_pv1_upper_bound_gr0]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv1_gr0
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv1_lower_bound_gr0]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv1_gr0
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_pv2_upper_bound_gr0]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv2_gr0
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv2_lower_bound_gr0]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv2_gr0
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_pv3_upper_bound_gr0]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv3_gr0
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv3_lower_bound_gr0]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv3_gr0
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_pv1_upper_bound_gr1]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv1_gr1
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv1_lower_bound_gr1]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv1_gr1
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_pv2_upper_bound_gr1]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv2_gr1
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv2_lower_bound_gr1]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv2_gr1
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_pv3_upper_bound_gr1]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv3_gr1
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv3_lower_bound_gr1]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv3_gr1
    bound_type = lower
    bound_value = -1
  [../]
[]


[Variables]
  # concentration
  [./c1]
    order = FIRST
    family = LAGRANGE
  [../]
  [./c2]
    order = FIRST
    family = LAGRANGE
  [../]

# phase concentration c1 in matrix
  [./c1m]
    order = FIRST
    family = LAGRANGE
  [../]
# phase concentration c1 in pv1
  [./c1pv1]
    order = FIRST
    family = LAGRANGE
  [../]
  # phase concentration c1 in pv2
  [./c1pv2]
    order = FIRST
    family = LAGRANGE
  [../]
 # phase concentration c1 in pv3
 [./c1pv3]
  order = FIRST
  family = LAGRANGE
[../]

# phase concentration c2 in matrix
  [./c2m]
    order = FIRST
    family = LAGRANGE
  [../]
  # phase concentration c2 in pv1
  [./c2pv1]
    order = FIRST
    family = LAGRANGE
  [../]
  # phase concentration c2 in pv2
  [./c2pv2]
    order = FIRST
    family = LAGRANGE
  [../]
 # phase concentration c2 in pv3
 [./c2pv3]
  order = FIRST
  family = LAGRANGE
[../]

  # order parameter m
  [./eta_m_gr0]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
  # order parameter pv1
  [./eta_pv1_gr0]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
  # order parameter pv2
  [./eta_pv2_gr0]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
# order parameter pv3
  [./eta_pv3_gr0]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
   # order parameter m
  [./eta_m_gr1]
    order = FIRST
    family = LAGRANGE
    block = 1
  [../]
  # order parameter pv1
  [./eta_pv1_gr1]
    order = FIRST
    family = LAGRANGE
    block = 1
  [../]
  # order parameter pv2
  [./eta_pv2_gr1]
    order = FIRST
    family = LAGRANGE
    block = 1
  [../]
# order parameter pv3
  [./eta_pv3_gr1]
    order = FIRST
    family = LAGRANGE
    block = 1
  [../]

  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Functions]
  [./bc_func]
    type = ParsedFunction
    expression = sin(alpha*pi*x)
    symbol_names = alpha
    symbol_values = 16
  [../]
  [./gb_scale_fn]
    type = ParsedFunction
    expression = '1 + (gb_factor - 1)*0.5*(tanh((w/2 - abs(x - x0))/delta) + 1)'
    symbol_names = 'x0          w     delta   gb_factor'
    symbol_values = '350.0     20.0      1       5.0'
  [../]
[]

[ICs]
  [./eta_pv1_gr0]
    variable = eta_pv1_gr0
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 324
    block = 0
  [../]
  [./eta_pv2_gr0]
    variable = eta_pv2_gr0
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 230	
    block = 0
  [../]
  [./eta_pv3_gr0]
    variable = eta_pv3_gr0
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 307	
    block = 0
  [../]
   [./eta_pv1_gr1]
    variable = eta_pv1_gr1
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 324
    block = 1
  [../]
  [./eta_pv2_gr1]
    variable = eta_pv2_gr1
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 230
    block = 1	
  [../]
  [./eta_pv3_gr1]
    variable = eta_pv3_gr1
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 307
    block = 1	
  [../]
  [./c1]
    variable = c1
    type = RandomIC
    min = 0.01775	
    max = 0.03025
    seed = 89	
  [../]
  [./c2]
    variable = c2
    type = RandomIC
    min = 0.03175	
    max = 0.04425
    seed = 89	
  [../]

[]


[Materials]
 # local free energy of the phase m in gr0
   [./fm_gr0]
    type = DerivativeParsedMaterial
    property_name = fc_m_gr0
    coupled_variables = 'c1m c2m'
    expression = '50.0*((c1m-0.0161)^2+2*(c2m-0.00723)^2)'
    block = 0
  [../]
  # Elastic energy of the phase m in gr0
  [./elastic_free_energy_m_gr0]
    type = ElasticEnergyMaterial
    base_name = phasem_gr0
    property_name = fe_m_gr0
    coupled_variables = ' '
    block = 0
  [../]
  # Total free energy of the phase m in gr0
  [./Total_energy_m_gr0]
    type = DerivativeSumMaterial
    property_name = Fm_gr0
    sum_materials = 'fc_m_gr0 fe_m_gr0'
    coupled_variables = 'c1m c2m'
    block = 0
  [../]
  # local free energy of the phase m in gr1
[./fm_gr1]
    type = DerivativeParsedMaterial
    property_name = fc_m_gr1
    coupled_variables = 'c1m c2m'
    expression = '50.0*((c1m-0.0161)^2+2*(c2m-0.00723)^2)'
    block = 1
  [../]
  # Elastic energy of the phase m in gr1
  [./elastic_free_energy_m_gr1]
    type = ElasticEnergyMaterial
    base_name = phasem_gr1
    property_name = fe_m_gr1
    coupled_variables = ' '
    block = 1
  [../]
  # Total free energy of the phase m in gr1
  [./Total_energy_m_gr1]
    type = DerivativeSumMaterial
    property_name = Fm_gr1
    sum_materials = 'fc_m_gr1 fe_m_gr1'
    coupled_variables = 'c1m c2m'
    block = 1
  [../]
  # local free energy of the phase pv1 in gr0
  [./fc_pv1_gr0]
    type = DerivativeParsedMaterial
    property_name = fc_pv1_gr0
    coupled_variables = 'c1pv1 c2pv1'
    expression = '50.0*((c1pv1-0.000727)^2+2*(c2pv1-0.196)^2)'
    block = 0
  [../]
  # Elastic energy of the phase pv1 in gr0
  [./elastic_free_energy_pv1_gr0]
    type = ElasticEnergyMaterial
    base_name = phasepv1_gr0
    property_name = fe_pv1_gr0
    coupled_variables = ' '
    block = 0
  [../]
  # Total free energy of the phase pv1 in gr0
  [./Total_energy_pv1_gr0]
    type = DerivativeSumMaterial
    property_name = Fpv1_gr0
    sum_materials = 'fc_pv1_gr0 fe_pv1_gr0'
    coupled_variables = 'c1pv1 c2pv1'
    block = 0
  [../]
  # local free energy of the phase pv1 in gr1
   [./fc_pv1_gr1]
    type = DerivativeParsedMaterial
    property_name = fc_pv1_gr1
    coupled_variables = 'c1pv1 c2pv1'
    expression = '50.0*((c1pv1-0.000727)^2+2*(c2pv1-0.196)^2)'
    block = 1
  [../]
  # Elastic energy of the phase pv1 in gr1
  [./elastic_free_energy_pv1_gr1]
    type = ElasticEnergyMaterial
    base_name = phasepv1_gr1
    property_name = fe_pv1_gr1
    coupled_variables = ' '
    block = 1
  [../]
  # Total free energy of the phase pv1 in gr1
  [./Total_energy_pv1_gr1]
    type = DerivativeSumMaterial
    property_name = Fpv1_gr1
    sum_materials = 'fc_pv1_gr1 fe_pv1_gr1'
    coupled_variables = 'c1pv1 c2pv1'
    block = 1
  [../]

   # local free energy of the phase pv2 in gr0
  [./f2_gr0]
    type = DerivativeParsedMaterial
    property_name = fc_pv2_gr0
    coupled_variables = 'c1pv2 c2pv2'
    expression = '50.0*((c1pv2-0.187)^2+2*(c2pv2-0.0157)^2)'
    block = 0
  [../]
  # Elastic energy of the phase pv2 in gr0
  [./elastic_free_energy_pv2_gr0]
    type = ElasticEnergyMaterial
    base_name = phasepv2_gr0
    property_name = fe_pv2_gr0
    coupled_variables = ' '
    block = 0
  [../]
  # Total free energy of the phase pv2 in gr0
  [./Total_energy_pv2_gr0]
    type = DerivativeSumMaterial
    property_name = Fpv2_gr0
    sum_materials = 'fc_pv2_gr0 fe_pv2_gr0'
    coupled_variables = 'c1pv2 c2pv2'
    block = 0
  [../]
 
  # local free energy of the phase pv2 in gr1
    [./f2_gr1]
    type = DerivativeParsedMaterial
    property_name = fc_pv2_gr1
    coupled_variables = 'c1pv2 c2pv2'
    expression = '50.0*((c1pv2-0.187)^2+2*(c2pv2-0.0157)^2)'
    block = 1
  [../]
  # Elastic energy of the phase pv2 in gr1
  [./elastic_free_energy_pv2_gr1]
    type = ElasticEnergyMaterial
    base_name = phasepv2_gr1
    property_name = fe_pv2_gr1
    coupled_variables = ' '
    block = 1
  [../]
  # Total free energy of the phase pv2 in gr1
  [./Total_energy_pv2_gr1]
    type = DerivativeSumMaterial
    property_name = Fpv2_gr1
    sum_materials = 'fc_pv2_gr1 fe_pv2_gr1'
    coupled_variables = 'c1pv2 c2pv2'
    block = 1
  [../]
 
  # local free energy of the phase pv3 in gr0
 [./f3_gr0]
    type = DerivativeParsedMaterial
    property_name = fc_pv3_gr0
    coupled_variables = 'c1pv3 c2pv3'
    expression = '50.0*((c1pv3-0.000727)^2+2*(c2pv3-0.196)^2)'
    block = 0
  [../]
    # Elastic energy of the phase pv3 in gr0
  [./elastic_free_energy_pv3_gr0]
    type = ElasticEnergyMaterial
    base_name = phasepv3_gr0
    property_name = fe_pv3_gr0
    coupled_variables = ' '
    block = 0
  [../]
    # Total free energy of the phase pv3 in gr0
  [./Total_energy_pv3_gr0]
    type = DerivativeSumMaterial
    property_name = Fpv3_gr0
    sum_materials = 'fc_pv3_gr0 fe_pv3_gr0'
    coupled_variables = 'c1pv3 c2pv3'
    block = 0
  [../]
  # local free energy of the phase pv3 in gr1
  [./f3_gr1]
    type = DerivativeParsedMaterial
    property_name = fc_pv3_gr1
    coupled_variables = 'c1pv3 c2pv3'
    expression = '50.0*((c1pv3-0.000727)^2+2*(c2pv3-0.196)^2)'
    block = 1
  [../]
    # Elastic energy of the phase pv3 in gr1
  [./elastic_free_energy_pv3_gr1]
    type = ElasticEnergyMaterial
    base_name = phasepv3_gr1
    property_name = fe_pv3_gr1
    coupled_variables = ' '
    block = 1
  [../]
    # Total free energy of the phase pv3 in gr1 
  [./Total_energy_pv3_gr1]
    type = DerivativeSumMaterial
    property_name = Fpv3_gr1
    sum_materials = 'fc_pv3_gr1 fe_pv3_gr1'
    coupled_variables = 'c1pv3 c2pv3'
    block = 1
  [../]
  # Switching functions for each phase
  # hm(eta_pv1, eta_pv2, eta_m)
  [./hm_gr0]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_m_gr0
    all_etas = 'eta_pv1_gr0 eta_pv2_gr0 eta_pv3_gr0 eta_m_gr0'
    h_name = hm_gr0
    block = 0
  [../]
  # hpv1(eta_pv1_gr0, eta_pv2_gr0,eta_pv3_gr0, eta_m_gr0)
  [./hpv1_gr0]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv1_gr0
    all_etas = 'eta_pv1_gr0 eta_pv2_gr0 eta_pv3_gr0 eta_m_gr0'
    h_name = hpv1_gr0
    block = 0
  [../]
  # hpv2(eta_pv1_gr0, eta_pv2_gr0,eta_pv3_gr0, eta_m_gr0)
  [./hpv2_gr0]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv2_gr0
    all_etas = 'eta_pv1_gr0 eta_pv2_gr0 eta_pv3_gr0 eta_m_gr0'
    h_name = hpv2_gr0
    block = 0
  [../]
  # hpv3(eta_pv1_gr0, eta_pv2_gr0, eta_pv3_gr0, eta_m_gr0)
[./hpv3_gr0]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv3_gr0
    all_etas = 'eta_pv1_gr0 eta_pv2_gr0 eta_pv3_gr0 eta_m_gr0'
    h_name = hpv3_gr0
    block = 0
[../]
  # hm(eta_pv1_gr1, eta_pv2_gr1, eta_pv3_gr1, eta_m_gr1)
[./hm_gr1]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_m_gr1
    all_etas = 'eta_pv1_gr1 eta_pv2_gr1 eta_pv3_gr1 eta_m_gr1'
    h_name = hm_gr1
    block = 1
  [../]
  # hpv1(eta_pv1_gr1, eta_pv2_gr1,eta_pv3_gr1, eta_m_gr1)
  [./hpv1_gr1]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv1_gr1
    all_etas = 'eta_pv1_gr1 eta_pv2_gr1 eta_pv3_gr1 eta_m_gr1'
    h_name = hpv1_gr1
    block = 1
  [../]
  # hpv2(eta_pv1_gr1, eta_pv2_gr1,eta_pv3_gr1, eta_m_gr1)
  [./hpv2_gr1]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv2_gr1
    all_etas = 'eta_pv1_gr1 eta_pv2_gr1 eta_pv3_gr1 eta_m_gr1'
    h_name = hpv2_gr1
    block = 1
  [../]
  # hpv3(eta_pv1_gr1, eta_pv2_gr1, eta_pv3_gr1, eta_m_gr1)
[./hpv3_gr1]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv3_gr1
    all_etas = 'eta_pv1_gr1 eta_pv2_gr1 eta_pv3_gr1 eta_m_gr1'
    h_name = hpv3_gr1
    block = 1
[../]
  
  # Coefficients for diffusion equation
  # --- self terms for c1: D11
  [./D11hm_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D11_gr0 hm_gr0'
    expression = gb_scale_aux*(D11_gr0*hm_gr0)
    property_name = D11hm_gr0
    block = 0
  [../]
  [./D11hpv1_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D11_gr0 hpv1_gr0'
    expression = gb_scale_aux*(D11_gr0*hpv1_gr0)
    property_name = D11hpv1_gr0
    block = 0
  [../]
  [./D11hpv2_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D11_gr0 hpv2_gr0'
    expression = gb_scale_aux*(D11_gr0*hpv2_gr0)
    property_name = D11hpv2_gr0
    block = 0
  [../]
  [./D11hpv3_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D11_gr0 hpv3_gr0'
    expression = gb_scale_aux*(D11_gr0*hpv3_gr0)
    property_name = D11hpv3_gr0
    block = 0
  [../]

# --- cross term for c1: D12 
  [./D12hm_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D12_gr0 hm_gr0'
    expression = gb_scale_aux*(D12_gr0*hm_gr0)  
    property_name = D12hm_gr0
    block = 0
  [../]
  [./D12hpv1_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D12_gr0 hpv1_gr0'
    expression = gb_scale_aux*(D12_gr0*hpv1_gr0)
    property_name = D12hpv1_gr0
    block = 0
  [../]
  [./D12hpv2_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D12_gr0 hpv2_gr0'
    expression = gb_scale_aux*(D12_gr0*hpv2_gr0)
    property_name = D12hpv2_gr0
    block = 0
  [../]
  [./D12hpv3_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D12_gr0 hpv3_gr0'
    expression = gb_scale_aux*(D12_gr0*hpv3_gr0)
    property_name = D12hpv3_gr0
    block = 0
  [../]

 #--- self terms for c2: D22
  [./D22hm_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D22_gr0 hm_gr0'
    expression = gb_scale_aux*(D22_gr0*hm_gr0)
    property_name = D22hm_gr0
    block = 0
  [../]
  [./D22hpv1_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D22_gr0 hpv1_gr0'
    expression = gb_scale_aux*(D22_gr0*hpv1_gr0)
    property_name = D22hpv1_gr0
    block = 0
  [../]
  [./D22hpv2_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D22_gr0 hpv2_gr0'
    expression = gb_scale_aux*(D22_gr0*hpv2_gr0)
    property_name = D22hpv2_gr0
    block = 0
  [../]
  [./D22hpv3_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D22_gr0 hpv3_gr0'
    expression = gb_scale_aux*(D22_gr0*hpv3_gr0)
    property_name = D22hpv3_gr0
    block = 0
  [../]
# --- cross term for c2: D21
  [./D21hm_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D21_gr0 hm_gr0'
    expression = gb_scale_aux*(D21_gr0*hm_gr0)
    property_name = D21hm_gr0
    block = 0
  [../]
  [./D21hpv1_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D21_gr0 hpv1_gr0'
    expression = gb_scale_aux*(D21_gr0*hpv1_gr0)
    property_name = D21hpv1_gr0
    block = 0
  [../]
  [./D21hpv2_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D21_gr0 hpv2_gr0'
    expression = gb_scale_aux*(D21_gr0*hpv2_gr0)
    property_name = D21hpv2_gr0
    block = 0
  [../]
  [./D21hpv3_gr0]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D21_gr0 hpv3_gr0'
    expression = gb_scale_aux*(D21_gr0*hpv3_gr0)
    property_name = D21hpv3_gr0
    block = 0
  [../]
# --- self terms for c1: D11
  [./D11hm_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D11_gr1 hm_gr1'
    expression = gb_scale_aux*(D11_gr1*hm_gr1)
    property_name = D11hm_gr1
    block = 1
  [../]
  [./D11hpv1_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D11_gr1 hpv1_gr1'
    expression = gb_scale_aux*(D11_gr1*hpv1_gr1)
    property_name = D11hpv1_gr1
    block = 1
  [../]
  [./D11hpv2_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D11_gr1 hpv2_gr1'
    expression = gb_scale_aux*(D11_gr1*hpv2_gr1)
    property_name = D11hpv2_gr1
    block = 1
  [../]
  [./D11hpv3_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D11_gr1 hpv3_gr1'
    expression = gb_scale_aux*(D11_gr1*hpv3_gr1)
    property_name = D11hpv3_gr1
    block = 1
  [../]
# --- cross term for c1: D12 
  [./D12hm_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D12_gr1 hm_gr1'
    expression = gb_scale_aux*(D12_gr1*hm_gr1)
    property_name = D12hm_gr1
    block = 1
  [../]
  [./D12hpv1_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D12_gr1 hpv1_gr1'
    expression = gb_scale_aux*(D12_gr1*hpv1_gr1)
    property_name = D12hpv1_gr1
    block = 1
  [../]
  [./D12hpv2_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D12_gr1 hpv2_gr1'
    expression = gb_scale_aux*(D12_gr1*hpv2_gr1)
    property_name = D12hpv2_gr1
    block = 1
  [../]
  [./D12hpv3_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D12_gr1 hpv3_gr1'
    expression = gb_scale_aux*(D12_gr1*hpv3_gr1)
    property_name = D12hpv3_gr1
    block = 1
  [../]
# --- self terms for c2: D22
  [./D22hm_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D22_gr1 hm_gr1'
    expression = gb_scale_aux*(D22_gr1*hm_gr1)
    property_name = D22hm_gr1
    block = 1
  [../]
  [./D22hpv1_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D22_gr1 hpv1_gr1'
    expression = gb_scale_aux*(D22_gr1*hpv1_gr1)
    property_name = D22hpv1_gr1
    block = 1
  [../]
  [./D22hpv2_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D22_gr1 hpv2_gr1'
    expression = gb_scale_aux*(D22_gr1*hpv2_gr1)
    property_name = D22hpv2_gr1
    block = 1
  [../]
  [./D22hpv3_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D22_gr1 hpv3_gr1'
    expression = gb_scale_aux*(D22_gr1*hpv3_gr1)
    property_name = D22hpv3_gr1
    block = 1
  [../]
# --- cross term for c2: D21
  [./D21hm_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D21_gr1 hm_gr1'
    expression = gb_scale_aux*(D21_gr1*hm_gr1)
    property_name = D21hm_gr1
    block = 1
  [../]
  [./D21hpv1_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D21_gr1 hpv1_gr1'
    expression = gb_scale_aux*(D21_gr1*hpv1_gr1)
    property_name = D21hpv1_gr1
    block = 1
  [../]
  [./D21hpv2_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D21_gr1 hpv2_gr1'
    expression = gb_scale_aux*(D21_gr1*hpv2_gr1)
    property_name = D21hpv2_gr1
    block = 1
  [../]
  [./D21hpv3_gr1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D21_gr1 hpv3_gr1'
    expression = gb_scale_aux*(D21_gr1*hpv3_gr1)
    property_name = D21hpv3_gr1
    block = 1
  [../]

# Barrier functions for each phase
  [./gm_gr0]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta_m_gr0
    function_name = gm_gr0
    block = 0
  [../]
  [./gpv1_gr0]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_pv1_gr0
    function_name = gpv1_gr0
    block = 0
  [../]
  [./gpv2_gr0]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_pv2_gr0
    function_name = gpv2_gr0
    block = 0
  [../]
  [./gpv3_gr0]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_pv3_gr0
    function_name = gpv3_gr0
    block = 0
  [../]
   [./gm_gr1]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta_m_gr1
    function_name = gm_gr1
    block = 1
  [../]
  [./gpv1_gr1]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_pv1_gr1
    function_name = gpv1_gr1
    block = 1
  [../]
  [./gpv2_gr1]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_pv2_gr1
    function_name = gpv2_gr1
    block = 1
  [../]
  [./gpv3_gr1]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_pv3_gr1
    function_name = gpv3_gr1
    block = 1
  [../]

 
  # constant properties
  [./constants_gr0]
    type = GenericConstantMaterial
    prop_names  = 'L_gr0    kappa_gr0  D11_gr0 D22_gr0 D12_gr0 D21_gr0 misfit_gr0     W_gr0'
    prop_values = '0.3       0.01        1        1     -0.5    -0.5      1            0.01'
    block = 0
  [../]
   [./constants_gr1]
    type = GenericConstantMaterial
    prop_names  = 'L_gr1    kappa_gr1  D11_gr1  D22_gr1  D12_gr1  D21_gr1  misfit_gr1     W_gr1'
    prop_values = '0.3         0.01      1        1        -0.5    -0.5         1         0.01'
    block = 1
  [../]

  #Mechanical properties
  [./Stiffness_phasem_gr0]
    type = ComputeElasticityTensor
    C_ijkl = '272.1 169 169 272.1 169 272.1 131 131 131' #Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of     
    base_name = phasem_gr0
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 0
  [../]
   [./Stiffness_phasem_gr1]
    type = ComputeElasticityTensor
    C_ijkl = '272.1 169 169 272.1 169 272.1 131 131 131' #Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of     
    base_name = phasem_gr1
    fill_method = symmetric9
     euler_angle_1 = 0
    euler_angle_2  = 0
    euler_angle_3  = 0
    block = 1
  [../]
  [./Stiffness_phasepv1_gr0]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'#Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of    
    base_name = phasepv1_gr0
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 0
  [../]
  [./Stiffness_phasepv1_gr1]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'#Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of    
    base_name = phasepv1_gr1
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 1
  [../]
  [./Stiffness_phasepv2_gr0]
    type = ComputeElasticityTensor
    C_ijkl = '243 154.8 154.8 243 154.8 243 132.3 132.3 132.3'
    base_name = phasepv2_gr0
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 0
  [../]
  [./Stiffness_phasepv2_gr1]
    type = ComputeElasticityTensor
    C_ijkl = '243 154.8 154.8 243 154.8 243 132.3 132.3 132.3'
    base_name = phasepv2_gr1
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 1
  [../]
  [./Stiffness_phasepv3_gr0]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
    base_name = phasepv3_gr0
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 0
  [../]
  [./Stiffness_phasepv3_gr1]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
    base_name = phasepv3_gr1
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 1
  [../]

  [./stress_phasepv1_gr0]
    type = ComputeLinearElasticStress
    base_name = phasepv1_gr0 
    block = 0
  [../]
  [./stress_phasepv1_gr1]
    type = ComputeLinearElasticStress
    base_name = phasepv1_gr1
    block = 1
  [../]
  [./stress_phasepv2_gr0]
    type = ComputeLinearElasticStress
    base_name = phasepv2_gr0
    block = 0
  [../]
  [./stress_phasepv2_gr1]
    type = ComputeLinearElasticStress
    base_name = phasepv2_gr1
    block = 1
  [../]
  [./stress_phasepv3_gr0]
    type = ComputeLinearElasticStress
    base_name = phasepv3_gr0
    block = 0
  [../]
  [./stress_phasepv3_gr1]
    type = ComputeLinearElasticStress
    base_name = phasepv3_gr1
    block = 1
  [../]
  [./stress_phasem_gr0]
    type = ComputeLinearElasticStress
    base_name = phasem_gr0
    block = 0
  [../]
  [./stress_phasem_gr1]
    type = ComputeLinearElasticStress
    base_name = phasem_gr1
    block = 1
  [../]

  [./strain_phasem_gr0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasem_gr0
    block = 0
  [../]
  [./strain_phasem_gr1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasem_gr1
    block = 1
  [../]
  [./strain_phasepv1_gr0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv1_gr0 
    eigenstrain_names = eigen_strainpv1_gr0
    block = 0
  [../]
  [./strain_phasepv1_gr1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv1_gr1
    eigenstrain_names = eigen_strainpv1_gr1
    block = 1
  [../]
  [./strain_phasepv2_gr0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv2_gr0
    eigenstrain_names = eigen_strainpv2_gr0
    block = 0
  [../]
  [./strain_phasepv2_gr1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv2_gr1
    eigenstrain_names = eigen_strainpv2_gr1
    block = 1
  [../]
  [./strain_phasepv3_gr0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv3_gr0
    eigenstrain_names = eigen_strainpv3_gr0
    block = 0
  [../]
  [./strain_phasepv3_gr1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv3_gr1
    eigenstrain_names = eigen_strainpv3_gr1
    block = 1
  [../]

  [./eigen_strainpv1_gr0]
    type = ComputeEigenstrain
    base_name = phasepv1_gr0
    eigen_base = '0.0067 0.0206 0 0 0 0'
    prefactor = misfit_gr0
    eigenstrain_name = eigen_strainpv1_gr0
    block = 0
  [../]
  [./eigen_strainpv1_gr1]
    type = ComputeEigenstrain
    base_name = phasepv1_gr1
    eigen_base = '0.0067 0.0206 0 0 0 0'
    prefactor = misfit_gr1
    eigenstrain_name = eigen_strainpv1_gr1
    block = 1
  [../]

  [./eigen_strainpv2_gr0]
    type = ComputeEigenstrain
    base_name = phasepv2_gr0
    eigen_base = '-0.003 -0.003 0 0 0 0'
    prefactor = misfit_gr0
    eigenstrain_name = eigen_strainpv2_gr0
    block = 0
  [../]
   [./eigen_strainpv2_gr1]
    type = ComputeEigenstrain
    base_name = phasepv2_gr1
    eigen_base = '-0.003 -0.003 0 0 0 0'
    prefactor = misfit_gr1
    eigenstrain_name = eigen_strainpv2_gr1
    block = 1
  [../]
  [./eigen_strainpv3_gr0]
    type = ComputeEigenstrain
    base_name = phasepv3_gr0
    eigen_base = '0.0206 0.0067 0 0 0 0'
    prefactor = misfit_gr0
    eigenstrain_name = eigen_strainpv3_gr0
    block = 0
  [../]
  [./eigen_strainpv3_gr1]
    type = ComputeEigenstrain
    base_name = phasepv3_gr1
    eigen_base =  '0.0206 0.0067 0 0 0 0'
    prefactor = misfit_gr1
    eigenstrain_name = eigen_strainpv3_gr1
    block = 1
  [../]


  # Generate the global stress from the phase stresses
  [./global_stress_gr0]
    type = MultiPhaseStressMaterial
    phase_base = 'phasepv1_gr0 phasepv2_gr0 phasepv3_gr0 phasem_gr0'
    h          = 'hpv1_gr0     hpv2_gr0   hpv3_gr0   hm_gr0'
    block = 0
  [../]
   [./global_stress_gr1]
    type = MultiPhaseStressMaterial
    phase_base = 'phasepv1_gr1 phasepv2_gr1 phasepv3_gr1 phasem_gr1'
    h          = 'hpv1_gr1     hpv2_gr1   hpv3_gr1   hm_gr1'
    block = 1
  [../]

  [./global_strain_gr0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    block = 0
  [../]
  [./global_strain_gr1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    block = 1
  [../]
[]

[Kernels] 
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
    incremental = true
  [../]
  # Kernels for Allen-Cahn equation for eta_pv1_gr0
  [./deta_pv1_dt_gr0]
    type = TimeDerivative
    variable = eta_pv1_gr0
    block = 0
  [../]
  [./ACBulkFpv1_gr0]
    type = KKSMultiACBulkF
    variable  = eta_pv1_gr0
    Fj_names  = 'Fpv1_gr0 Fpv2_gr0 Fpv3_gr0 Fm_gr0'
    hj_names  = 'hpv1_gr0 hpv2_gr0 hpv3_gr0 hm_gr0'
    gi_name   = gpv1_gr0
    eta_i     = eta_pv1_gr0
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m c2pv1 c2pv2 c2pv3 c2m eta_pv2_gr0 eta_pv3_gr0 eta_m_gr0'
    mob_name = L_gr0
    block = 0
  [../]
  [./ACBulkCpv1_c1_gr0]
    type = KKSMultiACBulkC
    variable  = eta_pv1_gr0
    Fj_names  = 'Fpv1_gr0 Fpv2_gr0 Fpv3_gr0 Fm_gr0'
    hj_names  = 'hpv1_gr0 hpv2_gr0 hpv3_gr0 hm_gr0'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1m'
    eta_i     = eta_pv1_gr0
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2m eta_pv2_gr0 eta_pv3_gr0 eta_m_gr0'
    mob_name = L_gr0
    block = 0
  [../]
  [./ACBulkCpv1_c2_gr0]
    type = KKSMultiACBulkC
    variable  = eta_pv1_gr0
    Fj_names  = 'Fpv1_gr0 Fpv2_gr0 Fpv3_gr0 Fm_gr0'
    hj_names  = 'hpv1_gr0 hpv2_gr0 hpv3_gr0 hm_gr0'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2m'
    eta_i     = eta_pv1_gr0
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m eta_pv2_gr0 eta_pv3_gr0 eta_m_gr0'
    mob_name = L_gr0
    block = 0
  [../]
  [./ACInterfacepv1_gr0]
    type = ACInterface
    variable = eta_pv1_gr0
    kappa_name = kappa_gr0
    mob_name = L_gr0
    block = 0
  [../]
    # Kernels for Allen-Cahn equation for eta_pv1_gr1
  [./deta_pv1_dt_gr1]
    type = TimeDerivative
    variable = eta_pv1_gr1
    block = 1
  [../]
  [./ACBulkFpv1_gr1]
    type = KKSMultiACBulkF
    variable  = eta_pv1_gr1
    Fj_names  = 'Fpv1_gr1 Fpv2_gr1 Fpv3_gr1 Fm_gr1'
    hj_names  = 'hpv1_gr1 hpv2_gr1 hpv3_gr1 hm_gr1'
    gi_name   = gpv1_gr1
    eta_i     = eta_pv1_gr1
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m c2pv1 c2pv2 c2pv3 c2m eta_pv2_gr1 eta_pv3_gr1 eta_m_gr1'
    mob_name = L_gr1
    block = 1
  [../]
  [./ACBulkCpv1_c1_gr1]
    type = KKSMultiACBulkC
    variable  = eta_pv1_gr1
    Fj_names  = 'Fpv1_gr1 Fpv2_gr1 Fpv3_gr1 Fm_gr1'
    hj_names  = 'hpv1_gr1 hpv2_gr1 hpv3_gr1 hm_gr1'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1m'
    eta_i     = eta_pv1_gr1
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2m eta_pv2_gr1 eta_pv3_gr1 eta_m_gr1'
    mob_name = L_gr1
    block = 1
  [../]
  [./ACBulkCpv1_c2_gr1]
    type = KKSMultiACBulkC
    variable  = eta_pv1_gr1
    Fj_names  = 'Fpv1_gr1 Fpv2_gr1 Fpv3_gr1 Fm_gr1'
    hj_names  = 'hpv1_gr1 hpv2_gr1 hpv3_gr1 hm_gr1'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2m'
    eta_i     = eta_pv1_gr1
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m eta_pv2_gr1 eta_pv3_gr1 eta_m_gr1'
    mob_name = L_gr1
    block = 1
  [../]
  [./ACInterfacepv1_gr1]
    type = ACInterface
    variable = eta_pv1_gr1
    kappa_name = kappa_gr1
    mob_name = L_gr1
    block = 1
  [../]

  # Kernels for Allen-Cahn equation for eta_pv2_gr0
  [./deta_pv2_dt_gr0]
    type = TimeDerivative
    variable = eta_pv2_gr0
    block = 0
  [../]
  [./ACBulkFpv2_gr0]
    type = KKSMultiACBulkF
    variable  = eta_pv2_gr0
    Fj_names  = 'Fpv1_gr0 Fpv2_gr0 Fpv3_gr0 Fm_gr0'
    hj_names  = 'hpv1_gr0 hpv2_gr0 hpv3_gr0 hm_gr0'
    gi_name   = gpv2_gr0
    eta_i     = eta_pv2_gr0
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m c2pv1 c2pv2 c2pv3 c2m eta_pv1_gr0 eta_pv3_gr0 eta_m_gr0'
    mob_name = L_gr0
    block = 0
  [../]
  [./ACBulkCpv2_c1_gr0]
    type = KKSMultiACBulkC
    variable  = eta_pv2_gr0
    Fj_names  = 'Fpv1_gr0 Fpv2_gr0 Fpv3_gr0 Fm_gr0'
    hj_names  = 'hpv1_gr0 hpv2_gr0 hpv3_gr0 hm_gr0'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1m'
    eta_i     = eta_pv2_gr0
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2m eta_pv1_gr0 eta_pv3_gr0 eta_m_gr0'
    mob_name = L_gr0
    block = 0
  [../]
  [./ACBulkCpv2_c2_gr0]
    type = KKSMultiACBulkC
    variable  = eta_pv2_gr0
    Fj_names  = 'Fpv1_gr0 Fpv2_gr0 Fpv3_gr0 Fm_gr0'
    hj_names  = 'hpv1_gr0 hpv2_gr0 hpv3_gr0 hm_gr0'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2m'
    eta_i     = eta_pv2_gr0
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m eta_pv1_gr0 eta_pv3_gr0 eta_m_gr0'
    mob_name = L_gr0
    block = 0
  [../]
  [./ACInterfacepv2_gr0]
    type = ACInterface
    variable = eta_pv2_gr0
    kappa_name = kappa_gr0
    mob_name = L_gr0
    block = 0
  [../]

 # Kernels for Allen-Cahn equation for eta_pv2_gr1
  [./deta_pv2_dt_gr1]
    type = TimeDerivative
    variable = eta_pv2_gr1
    block = 1
  [../]
  [./ACBulkFpv2_gr1]
    type = KKSMultiACBulkF
    variable  = eta_pv2_gr1
    Fj_names  = 'Fpv1_gr1 Fpv2_gr1 Fpv3_gr1 Fm_gr1'
    hj_names  = 'hpv1_gr1 hpv2_gr1 hpv3_gr1 hm_gr1'
    gi_name   = gpv2_gr1
    eta_i     = eta_pv2_gr1
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m c2pv1 c2pv2 c2pv3 c2m eta_pv1_gr1 eta_pv3_gr1 eta_m_gr1'
    mob_name = L_gr1
    block = 1
  [../]
  [./ACBulkCpv2_c1_gr1]
    type = KKSMultiACBulkC
    variable  = eta_pv2_gr1
    Fj_names  = 'Fpv1_gr1 Fpv2_gr1 Fpv3_gr1 Fm_gr1'
    hj_names  = 'hpv1_gr1 hpv2_gr1 hpv3_gr1 hm_gr1'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1m'
    eta_i     = eta_pv2_gr1
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2m eta_pv1_gr1 eta_pv3_gr1 eta_m_gr1'
    mob_name = L_gr1
    block = 1
  [../]
  [./ACBulkCpv2_c2_gr1]
    type = KKSMultiACBulkC
    variable  = eta_pv2_gr1
    Fj_names  = 'Fpv1_gr1 Fpv2_gr1 Fpv3_gr1 Fm_gr1'
    hj_names  = 'hpv1_gr1 hpv2_gr1 hpv3_gr1 hm_gr1'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2m'
    eta_i     = eta_pv2_gr1
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m eta_pv1_gr1 eta_pv3_gr1 eta_m_gr1'
    mob_name = L_gr1
    block = 1
  [../]
  [./ACInterfacepv2_gr1]
    type = ACInterface
    variable = eta_pv2_gr1
    kappa_name = kappa_gr1
    mob_name = L_gr1
    block = 1
  [../]

  # Kernels for Allen-Cahn equation for eta_pv3_gr0
  [./deta_pv3_dt_gr0]
    type = TimeDerivative
    variable = eta_pv3_gr0
    block = 0
  [../]
  [./ACBulkFpv3_gr0]
    type = KKSMultiACBulkF
    variable  = eta_pv3_gr0
    Fj_names  = 'Fpv1_gr0 Fpv2_gr0 Fpv3_gr0 Fm_gr0'
    hj_names  = 'hpv1_gr0 hpv2_gr0 hpv3_gr0 hm_gr0'
    gi_name   = gpv3_gr0
    eta_i     = eta_pv3_gr0
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m c2pv1 c2pv2 c2pv3 c2m eta_pv1_gr0 eta_pv2_gr0 eta_m_gr0'
    mob_name = L_gr0
    block = 0
  [../]
  [./ACBulkCpv3_c1_gr0]
    type = KKSMultiACBulkC
    variable  = eta_pv3_gr0
    Fj_names  = 'Fpv1_gr0 Fpv2_gr0 Fpv3_gr0 Fm_gr0'
    hj_names  = 'hpv1_gr0 hpv2_gr0 hpv3_gr0 hm_gr0'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1m'
    eta_i     = eta_pv3_gr0
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2m eta_pv1_gr0 eta_pv2_gr0 eta_m_gr0'
    mob_name = L_gr0
    block = 0
  [../]
  [./ACBulkCpv3_c2_gr0]
    type = KKSMultiACBulkC
    variable  = eta_pv3_gr0
    Fj_names  = 'Fpv1_gr0 Fpv2_gr0 Fpv3_gr0 Fm_gr0'
    hj_names  = 'hpv1_gr0 hpv2_gr0 hpv3_gr0 hm_gr0'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2m'
    eta_i     = eta_pv3_gr0
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m eta_pv1_gr0 eta_pv2_gr0 eta_m_gr0'
    mob_name = L_gr0
    block = 0
  [../]
  [./ACInterfacepv3_gr0]
    type = ACInterface
    variable = eta_pv3_gr0
    kappa_name = kappa_gr0
    mob_name = L_gr0
    block = 0
  [../]

 # Kernels for Allen-Cahn equation for eta_pv3_gr1
  [./deta_pv3_dt_gr1]
    type = TimeDerivative
    variable = eta_pv3_gr1
    block = 1
  [../]
  [./ACBulkFpv3_gr1]
    type = KKSMultiACBulkF
    variable  = eta_pv3_gr1
    Fj_names  = 'Fpv1_gr1 Fpv2_gr1 Fpv3_gr1 Fm_gr1'
    hj_names  = 'hpv1_gr1 hpv2_gr1 hpv3_gr1 hm_gr1'
    gi_name   = gpv3_gr1
    eta_i     = eta_pv3_gr1
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m c2pv1 c2pv2 c2pv3 c2m eta_pv1_gr1 eta_pv2_gr1 eta_m_gr1'
    mob_name = L_gr1
    block = 1
  [../]
  [./ACBulkCpv3_c1_gr1]
    type = KKSMultiACBulkC
    variable  = eta_pv3_gr1
    Fj_names  = 'Fpv1_gr1 Fpv2_gr1 Fpv3_gr1 Fm_gr1'
    hj_names  = 'hpv1_gr1 hpv2_gr1 hpv3_gr1 hm_gr1'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1m'
    eta_i     = eta_pv3_gr1
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2m eta_pv1_gr1 eta_pv2_gr1 eta_m_gr1'
    mob_name = L_gr1
    block = 1
  [../]
  [./ACBulkCpv3_c2_gr1]
    type = KKSMultiACBulkC
    variable  = eta_pv3_gr1
    Fj_names  = 'Fpv1_gr1 Fpv2_gr1 Fpv3_gr1 Fm_gr1'
    hj_names  = 'hpv1_gr1 hpv2_gr1 hpv3_gr1 hm_gr1'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2m'
    eta_i     = eta_pv3_gr1
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m eta_pv1_gr1 eta_pv2_gr1 eta_m_gr1'
    mob_name = L_gr1
    block = 1
  [../]
  [./ACInterfacepv3_gr1]
    type = ACInterface
    variable = eta_pv3_gr1
    kappa_name = kappa_gr1
    mob_name = L_gr1
    block = 1
  [../]

# Kernels for constraint equation |eta_pv1| + |eta_pv2| + eta_m = 1
  # eta3 is the nonlinear variable for the constraint equation
  [./eta_mreaction_gr0]
    type = MatReaction
    variable = eta_m_gr0
    reaction_rate = 1
    block = 0
 [../]
  [./eta_pv1reaction_gr0]
    type = MatReaction_abscouple
    variable = eta_m_gr0
    v = eta_pv1_gr0
    reaction_rate = 1
    block = 0
  [../]
  [./eta_pv2reaction_gr0]
    type = MatReaction_abscouple
    variable = eta_m_gr0
    v = eta_pv2_gr0
    reaction_rate = 1
    block = 0
  [../]
  [./eta_pv3reaction_gr0]
    type = MatReaction_abscouple
    variable = eta_m_gr0
    v = eta_pv3_gr0
    reaction_rate = 1
    block = 0
    [../]
   [./eta_mreaction_gr1]
    type = MatReaction
    variable = eta_m_gr1
    reaction_rate = 1
    block = 1
 [../]
  [./eta_pv1reaction_gr1]
    type = MatReaction_abscouple
    variable = eta_m_gr1
    v = eta_pv1_gr1
    reaction_rate = 1
    block = 1
  [../]
  [./eta_pv2reaction_gr1]
    type = MatReaction_abscouple
    variable = eta_m_gr1
    v = eta_pv2_gr1
    reaction_rate = 1
    block = 1
  [../]
  [./eta_pv3reaction_gr1]
    type = MatReaction_abscouple
    variable = eta_m_gr1
    v = eta_pv3_gr1
    reaction_rate = 1
    block = 1
    [../]
  [./one_gr0]
    type = BodyForce
    variable = eta_m_gr0
    value = -1.0
  [../]
  [./one_gr1]
    type = BodyForce
    variable = eta_m_gr1
    value = -1.0
  [../]

  #Kernels for diffusion equation of c1
  [./diff_time_c1]
    type = TimeDerivative
    variable = c1
  [../]
  [./diff_time_c2]
    type = TimeDerivative
    variable = c2
  [../]
   # ===== Grain 0 (block = 0) =====
  # ==== c1 equation: self terms (D11 ∇c1) ====
  [./diff_c1m_gr0]
    type = MatDiffusion
    variable = c1
    diffusivity = D11hm_gr0
    v = c1m
    block = 0
  [../]
  [./diff_c1pv1_gr0]
    type = MatDiffusion
    variable = c1
    diffusivity = D11hpv1_gr0
    v = c1pv1
    block = 0
  [../]
  [./diff_c1pv2_gr0]
    type = MatDiffusion
    variable = c1
    diffusivity = D11hpv2_gr0
    v = c1pv2
    block = 0
  [../]
  [./diff_c1pv3_gr0]
    type = MatDiffusion
    variable = c1
    diffusivity = D11hpv3_gr0
    v = c1pv3
    block = 0
  [../]

   # ===== Grain 0 (block = 0) =====
 # ==== c1 equation: cross terms (D12 ∇c2) ====
 [./diff_c1m_gr0_cross]
   type = MatDiffusion
    variable = c1
    diffusivity = D12hm_gr0
    v = c2m
    block = 0
  [../]
  [./diff_c1pv1_gr0_cross]
    type = MatDiffusion
    variable = c1
    diffusivity = D12hpv1_gr0
    v = c2pv1
    block = 0
  [../]
  [./diff_c1pv2_gr0_cross]
    type = MatDiffusion
    variable = c1
    diffusivity = D12hpv2_gr0
    v = c2pv2
    block = 0
  [../]
  [./diff_c1pv3_gr0_cross]
    type = MatDiffusion
    variable = c1
    diffusivity = D12hpv3_gr0
    v = c2pv3
    block = 0
  [../]

 # ===== Grain 0 (block = 0) =====
 # ==== c2 equation: cross terms (D22 ∇c2) ====
 [./diff_c2m_gr0]
    type = MatDiffusion
    variable = c2
    diffusivity = D22hm_gr0
    v = c2m
    block = 0
  [../]
  [./diff_c2pv1_gr0]
    type = MatDiffusion
    variable = c2
    diffusivity = D22hpv1_gr0
    v = c2pv1
    block = 0
  [../]
  [./diff_c2pv2_gr0]
    type = MatDiffusion
    variable = c2
    diffusivity = D22hpv2_gr0
    v = c2pv2
     block = 0
  [../]
  [./diff_c2pv3_gr0]
    type = MatDiffusion
    variable = c2
    diffusivity = D22hpv3_gr0
    v = c2pv3
    block = 0
  [../]
  
  # ===== Grain 0 (block = 0) =====
 # ==== c2 equation: cross terms (D21 ∇c1) ====
  [./diff_c2m_gr0_cross]
    type = MatDiffusion
    variable = c2
    diffusivity = D21hm_gr0
    v = c1m
    block = 0
  [../]
  [./diff_c2pv1_gr0_cross]
    type = MatDiffusion
    variable = c2
    diffusivity = D21hpv1_gr0
    v = c1pv1
    block = 0
  [../]
  [./diff_c2pv2_gr0_cross]
    type = MatDiffusion  
    variable = c2
    diffusivity = D21hpv2_gr0
    v = c1pv2
     block = 0
  [../]
  [./diff_c2pv3_gr0_cross]
    type = MatDiffusion
    variable = c2
    diffusivity = D21hpv3_gr0
    v = c1pv3
    block = 0
  [../]


   # ===== Grain 1 (block = 1) =====
  # ==== c1 equation: self terms (D11 ∇c1) ====
  [./diff_c1m_gr1]
    type = MatDiffusion
    variable = c1
    diffusivity = D11hm_gr1
    v = c1m
    block = 1
  [../] 
  [./diff_c1pv1_gr1]
    type = MatDiffusion
    variable = c1
    diffusivity = D11hpv1_gr1
    v = c1pv1
    block = 1
  [../]
  [./diff_c1pv2_gr1]
    type = MatDiffusion
    variable = c1
    diffusivity = D11hpv2_gr1
    v = c1pv2
     block = 1
  [../]
  [./diff_c1pv3_gr1]
    type = MatDiffusion
    variable = c1
    diffusivity = D11hpv3_gr1
    v = c1pv3
    block = 1
  [../]

   # ===== Grain 1 (block = 1) =====
 # ==== c1 equation: cross terms (D12 ∇c2) ====
 [./diff_c1m_gr1_cross]
    type = MatDiffusion
    variable = c1
    diffusivity = D12hm_gr1
    v = c2m
    block = 1 
  [../]
  [./diff_c1pv1_gr1_cross]
    type = MatDiffusion
    variable = c1
    diffusivity = D12hpv1_gr1
    v = c2pv1
    block = 1
  [../]
  [./diff_c1pv2_gr1_cross]
    type = MatDiffusion
    variable = c1
    diffusivity = D12hpv2_gr1
    v = c2pv2
     block = 1
  [../]
  [./diff_c1pv3_gr1_cross]
    type = MatDiffusion
    variable = c1
    diffusivity = D12hpv3_gr1
    v = c2pv3
    block = 1
  [../]
 
  # ===== Grain 1 (block = 1) =====
 # ==== c2 equation: self terms (D22 ∇c2) ====  
  [./diff_c2m_gr1]
    type = MatDiffusion
    variable = c2
    diffusivity = D22hm_gr1
    v = c2m
    block = 1
  [../]
  [./diff_c2pv1_gr1]
    type = MatDiffusion
    variable = c2
    diffusivity = D22hpv1_gr1
    v = c2pv1
    block = 1
  [../]
  [./diff_c2pv2_gr1]
    type = MatDiffusion
    variable = c2
    diffusivity = D22hpv2_gr1
    v = c2pv2
     block = 1
  [../]
  [./diff_c2pv3_gr1]
    type = MatDiffusion
    variable = c2
    diffusivity = D22hpv3_gr1
    v = c2pv3
    block = 1
  [../]
  # ===== Grain 1 (block = 1) =====
 # ==== c2 equation: cross terms (D21 ∇c1) ====
  [./diff_c2m_gr1_cross]
    type = MatDiffusion
    variable = c2
    diffusivity = D21hm_gr1
    v = c1m
    block = 1
  [../]
  [./diff_c2pv1_gr1_cross]
    type = MatDiffusion
    variable = c2
    diffusivity = D21hpv1_gr1
    v = c1pv1
    block = 1
  [../]
  [./diff_c2pv2_gr1_cross]
    type = MatDiffusion
    variable = c2
    diffusivity = D21hpv2_gr1
    v = c1pv2
     block = 1
  [../]
  [./diff_c2pv3_gr1_cross]
    type = MatDiffusion
    variable = c2
    diffusivity = D21hpv3_gr1
    v = c1pv3
    block = 1
  [../]
    
  # Phase concentration constraints
    [./chempot1m_pv1_gr0]
    type = KKSPhaseChemicalPotential
    variable = c1m
    cb       = c1pv1
    fa_name  = Fm_gr0
    fb_name  = Fpv1_gr0
    args_a   = c2m
    args_b   = c2pv1
    block = 0
  [../]
 [./chempot1m_pv2_gr0]
    type = KKSPhaseChemicalPotential
    variable = c1pv1
    cb       = c1pv2
    fa_name  = Fpv1_gr0
    fb_name  = Fpv2_gr0
    args_a   = c2pv1
    args_b   = c2pv2
    block = 0
  [../]
  [./chempot1m_pv3_gr0]
    type = KKSPhaseChemicalPotential
    variable = c1pv2
    cb       = c1pv3
    fa_name  = Fpv2_gr0
    fb_name  = Fpv3_gr0
    args_a   = c2pv2
    args_b   = c2pv3
    block = 0
  [../]
  [./chempot2m_pv1_gr0]
    type = KKSPhaseChemicalPotential
    variable = c2m
    cb       = c2pv1
    fa_name  = Fm_gr0
    fb_name  = Fpv1_gr0
    args_a   = c1m
    args_b   = c1pv1
    block = 0
  [../]
  [./chempot2m_pv2_gr0]
    type = KKSPhaseChemicalPotential
    variable = c2pv1
    cb       = c2pv2
    fa_name  = Fpv1_gr0
    fb_name  = Fpv2_gr0
    args_a   = c1pv1
    args_b   = c1pv2
    block = 0
  [../]
  [./chempot2m_pv3_gr0]
    type = KKSPhaseChemicalPotential
    variable = c2pv2
    cb       = c2pv3
    fa_name  = Fpv2_gr0
    fb_name  = Fpv3_gr0
    args_a   = c1pv2
    args_b   = c1pv3
     block = 0
  [../]
  

   [./chempot1m_pv1_gr1]
    type = KKSPhaseChemicalPotential
    variable = c1m
    cb       = c1pv1
    fa_name  = Fm_gr1
    fb_name  = Fpv1_gr1
    args_a   = c2m
    args_b   = c2pv1
    block = 1
  [../]
 [./chempot1m_pv2_gr1]
    type = KKSPhaseChemicalPotential
    variable = c1pv1
    cb       = c1pv2
    fa_name  = Fpv1_gr1
    fb_name  = Fpv2_gr1
    args_a   = c2pv1
    args_b   = c2pv2
    block = 1
  [../]
  [./chempot1m_pv3_gr1]
    type = KKSPhaseChemicalPotential
    variable = c1pv2
    cb       = c1pv3
    fa_name  = Fpv2_gr1
    fb_name  = Fpv3_gr1
    args_a   = c2pv2
    args_b   = c2pv3
    block = 1
  [../]
  [./chempot2m_pv1_gr1]
    type = KKSPhaseChemicalPotential
    variable = c2m
    cb       = c2pv1
    fa_name  = Fm_gr1
    fb_name  = Fpv1_gr1
    args_a   = c1m
    args_b   = c1pv1
    block = 1
  [../]
  [./chempot2m_pv2_gr1]
    type = KKSPhaseChemicalPotential
    variable = c2pv1
    cb       = c2pv2
    fa_name  = Fpv1_gr1
    fb_name  = Fpv2_gr1
    args_a   = c1pv1
    args_b   = c1pv2
    block = 1
  [../]
  [./chempot2m_pv3_gr1]
    type = KKSPhaseChemicalPotential
    variable = c2pv2
    cb       = c2pv3
    fa_name  = Fpv2_gr1
    fb_name  = Fpv3_gr1
    args_a   = c1pv2
    args_b   = c1pv3
    block = 1
  [../]
    
  [./phaseconcentration_c1pv1_gr0]
    type = KKSMultiPhaseConcentration
    variable = c1pv3
    cj = 'c1m c1pv1 c1pv2 c1pv3'
    hj_names = 'hm_gr0 hpv1_gr0 hpv2_gr0 hpv3_gr0'
    etas = 'eta_m_gr0 eta_pv1_gr0 eta_pv2_gr0 eta_pv3_gr0'
    c = c1
    block = 0
  [../]
  [./phaseconcentration_c2pv1_gr0]
    type = KKSMultiPhaseConcentration
    variable = c2pv3
    cj = 'c2m c2pv1 c2pv2 c2pv3'
    hj_names = 'hm_gr0 hpv1_gr0 hpv2_gr0 hpv3_gr0'
    etas = 'eta_m_gr0 eta_pv1_gr0 eta_pv2_gr0 eta_pv3_gr0'
    c = c2
    block = 0
  [../]
  [./phaseconcentration_c1pv1_gr1]
    type = KKSMultiPhaseConcentration
    variable = c1pv3
    cj = 'c1m c1pv1 c1pv2 c1pv3'
    hj_names = 'hm_gr1 hpv1_gr1 hpv2_gr1 hpv3_gr1'
    etas = 'eta_m_gr1 eta_pv1_gr1 eta_pv2_gr1 eta_pv3_gr1'
    c = c1
    block  = 1
  [../]
  [./phaseconcentration_c2pv1_gr1]
    type = KKSMultiPhaseConcentration
    variable = c2pv3
    cj = 'c2m c2pv1 c2pv2 c2pv3'
    hj_names = 'hm_gr1 hpv1_gr1 hpv2_gr1 hpv3_gr1'
    etas = 'eta_m_gr1 eta_pv1_gr1 eta_pv2_gr1 eta_pv3_gr1'
    c = c2
    block = 1
  [../]
  

[]

[AuxKernels]
  [./gb_scale_eval]
    type     = FunctionAux
    variable = gb_scale_aux
    function = gb_scale_fn
  [../]
  [temperature_gr0]
    type = FunctionAux
    variable = temperature_gr0
    function = '1023'
    execute_on = timestep_begin
    block = 0
  []
   [temperature_gr1]
    type = FunctionAux
    variable = temperature_gr1
    function = '1023'
    execute_on = timestep_begin
    block = 1
  []
  [./Energy_total]
    type = KKSMultiFreeEnergy
    Fj_names = 'Fpv1_gr0 Fpv2_gr0 Fm_gr0'
    hj_names = 'hpv1_gr0 hpv2_gr0 hm_gr0'
    gj_names = 'gpv1_gr0 gpv2_gr0 gm_gr0'
    variable = Energy
    w = 1
    interfacial_vars =  'eta_pv1_gr0  eta_pv2_gr0  eta_m_gr0'
    kappa_names =       'kappa_gr0 kappa_gr0 kappa_gr0'
    block = 0
  [../]
  [./stress_xx_gr0]
    type = RankTwoAux
    variable = stress_xx_gr0
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./stress_xx_gr1]
    type = RankTwoAux
    variable = stress_xx_gr1
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu            mumps            vinewtonrsls'

  l_max_its = 50
  nl_max_its = 25
  l_tol = 1.0e-3
  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-8

  end_time = 50000

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 5e-4
    cutback_factor = 0.75
    growth_factor = 1.2
    optimal_iterations = 20
  [../]

[]

[Preconditioning]
  active = 'full'
  [./full]
    type = SMP
    full = true
  [../]
  
[]

[Postprocessors]
   [./dofs]
     type = NumDOFs
   [../]
   [./h1_error]
     type = ElementH1Error
     variable = eta_pv1_gr0
     function = bc_func
     block = 0
   [../]
  [./pv1_area_gr0]
    type = ElementIntegralVariablePostprocessor_new2
    variable = eta_pv1_gr0
    block = 0
  [../]
   [./pv2_area_gr0]
    type = ElementIntegralVariablePostprocessor_new2
    variable = eta_pv2_gr0
    block = 0
  [../]
  [./pv3_area_gr0]
    type = ElementIntegralVariablePostprocessor_new2
    variable = eta_pv3_gr0
    block = 0
  [../]
   [./pv1_area_gr1]
    type = ElementIntegralVariablePostprocessor_new2
    variable = eta_pv1_gr1
    block = 1
  [../]
   [./pv2_area_gr1]
    type = ElementIntegralVariablePostprocessor_new2
    variable = eta_pv2_gr1
    block = 1
  [../]
  [./pv3_area_gr1]
    type = ElementIntegralVariablePostprocessor_new2
    variable = eta_pv3_gr1
    block = 1
  [../]
[]

[Outputs]
  exodus = true
  csv = true
  [./table]
    type = CSV
    execute_on = timestep_end
  [../]
  [./checkpoint]
     type = Checkpoint
     num_files = 10
     interval = 10
  [../]
[]
