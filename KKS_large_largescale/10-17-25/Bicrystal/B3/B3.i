# This test is for the multicomponent In718 alloy

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
   # Smooth GB band weight (dimensionless)
  [./gb_mask]
    family = MONOMIAL
    order  = CONSTANT
  [../]

  # Weighted fields (for integrals/averages)
  [./c2_gb_w]
    family = MONOMIAL
    order  = CONSTANT
  [../]
  [./dc2_gb_w]
    family = MONOMIAL
    order  = CONSTANT
  [../]
  [./h_g_pv1] 
  family=MONOMIAL 
  order=CONSTANT 
  [../]
  [./h_g_pv2] 
  family=MONOMIAL 
  order=CONSTANT 
  [../]
  [./h_g_pv3] 
  family=MONOMIAL 
  order=CONSTANT 
  [../]
  [./h_g_pv4] 
  family=MONOMIAL 
  order=CONSTANT 
  [../]
  [./g_aux]  
  family=MONOMIAL 
  order=CONSTANT 
  [../]
  [./one]    
  family=MONOMIAL 
  order=CONSTANT 
  [../]
   [./gb_scale_aux]
    family = MONOMIAL
    order  = FIRST
  [../]
  [./h_pv1_aux]
     family = MONOMIAL 
     order = CONSTANT 
  [../]
  [./h_pv2_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./h_pv3_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./h_pv4_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./h_m_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
   [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_h_pv1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_h_pv2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_h_pv3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_h_pv4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_h_m]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./Energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]
  [temperature]
    order = FIRST
    family = LAGRANGE
  []
[]

[Bounds]
  [./eta_pv1_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv1
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv1_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv1
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_pv2_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv2
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv2_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv2
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_pv3_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv3
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv3_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv3
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
  # phase concentration c1 in pv4
  [./c1pv4]
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
  # phase concentration c2 in pv4
  [./c2pv4]
    order = FIRST
    family = LAGRANGE
  [../]

  # order parameter m
  [./eta_m]
    order = FIRST
    family = LAGRANGE
  [../]
  # order parameter pv1
  [./eta_pv1]
    order = FIRST
    family = LAGRANGE
  [../]
  # order parameter pv2
  [./eta_pv2]
    order = FIRST
    family = LAGRANGE
  [../]
# order parameter pv3
[./eta_pv3]
  order = FIRST
  family = LAGRANGE
[../]
  # order parameter pv4
  [./eta_pv4]
    order = FIRST
    family = LAGRANGE
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
    symbol_values = '350.0     30.0     1     20'
  [../]

   # Smooth switch: ~0 before 3 h, ~1 after 3 h
  # Convert simulation time -> hours with t/6630.57
  #[./gate_3h]
   # type  = ParsedFunction
   # value = '0.5*(1 + tanh((t/6630.57 - 3.0)/1.0))'
  #[../]

  # (Optional) smooth tiny δ seed at the GB if you want to guarantee nucleation there
  [./eta4_tiny_gb_fn]
    type  = ParsedFunction
    value = '0.02*exp(-pow((x - 350)/1.0, 2))'            # ~5 nm half-width, amplitude 0.02
  [../] 
# ---- Nb (c2) enrichment at the grain boundary (x = 350 sim units) ----
  # Baseline values from Wang: x_Al ≈ 0.024, x_Nb ≈ 0.038 (adjust if your alloy differs).
  # dC_nb is the peak enrichment (keep modest to preserve mass balance).
  [./c2_gb_enrich_fn]
    type  = ParsedFunction
    value = '0.038 + 0.006*exp(-pow((x - 350)/5.0, 2))'   # center at GB, ~25 nm half-width
  [../]
[]

[ICs]
  [./eta_pv1]
    variable = eta_pv1
    type = RandomIC
    min = 0
    max = 0.1
    seed = 324
  [../]
  [./eta_pv2]
    variable = eta_pv2
    type = RandomIC
    min = 0
    max = 0.1
    seed = 230	
  [../]
  [./eta_pv3]
    variable = eta_pv3
    type = RandomIC
    min = 0
    max = 0.1
    seed = 307	
  [../]
 
  [./eta_pv4_zero]
    type=FunctionIC
    variable = eta_pv4
    function = eta4_tiny_gb_fn
  [../]

    [./c1_global_ic]
    type     = ConstantIC
    variable = c1
    value    = 0.024            # adjust to your alloy if different
  [../]
  [./c2_global_ic]
    type     = FunctionIC
    variable = c2
    function = c2_gb_enrich_fn
  [../]
 # ---------- KKS per-phase compositions ----------
  # Initialize all per-phase c2 fields from the same GB-enriched profile
  [./c2m_ic]
    type=FunctionIC
    variable=c2m
    function=c2_gb_enrich_fn
  [../]
  [./c2pv1_ic]
    type=FunctionIC
    variable=c2pv1
    function=c2_gb_enrich_fn
  [../]
  [./c2pv2_ic]
    type=FunctionIC
    variable=c2pv2
    function=c2_gb_enrich_fn
  [../]
  [./c2pv3_ic]
    type=FunctionIC
    variable=c2pv3
    function=c2_gb_enrich_fn
  [../]
  [./c2pv4_ic]
    type=FunctionIC
    variable=c2pv4
    function=c2_gb_enrich_fn
  [../]

  # Al-like per-phase (c1*): match the global baseline for a consistent start
  [./c1m_ic]
    type=ConstantIC
    variable=c1m
    value=0.024
  [../]
  [./c1pv1_ic]
    type=ConstantIC
    variable=c1pv1
    value=0.024
  [../]
  [./c1pv2_ic]
    type=ConstantIC
    variable=c1pv2
    value=0.024
  [../]
  [./c1pv3_ic]
    type=ConstantIC
    variable=c1pv3
    value=0.024
  [../]
  [./c1pv4_ic]
    type=ConstantIC
    variable=c1pv4
    value=0.024
  [../]
[]


[Materials]
  [./fm]
    type = DerivativeParsedMaterial
    property_name = fc_m
    coupled_variables = 'c1m c2m'
    expression = '50.0*((c1m-0.0161)^2+2*(c2m-0.00723)^2)'
  [../]
  # Elastic energy of the phase 0
  [./elastic_free_energy_m]
    type = ElasticEnergyMaterial
    base_name = phasem
    property_name = fe_m
    coupled_variables = ' '
  [../]
  # Total free energy of the phase 0
  [./Total_energy_m]
    type = DerivativeSumMaterial
    property_name = Fm
    sum_materials = 'fc_m fe_m'
    coupled_variables = 'c1m c2m'
  [../]

  [./fc_pv1]
    type = DerivativeParsedMaterial
    property_name = fc_pv1
    coupled_variables = 'c1pv1 c2pv1'
    expression = '50.0*((c1pv1-0.000727)^2+2*(c2pv1-0.196)^2)'
  [../]
  # Elastic energy of the phase 1
  [./elastic_free_energy_pv1]
    type = ElasticEnergyMaterial
    base_name = phasepv1
    property_name = fe_pv1
    coupled_variables = ' '
  [../]
  # Total free energy of the phase 1
  [./Total_energy_pv1]
    type = DerivativeSumMaterial
    property_name = Fpv1
    sum_materials = 'fc_pv1 fe_pv1'
    coupled_variables = 'c1pv1 c2pv1'
  [../]

  [./f2]
    type = DerivativeParsedMaterial
    property_name = fc_pv2
    coupled_variables = 'c1pv2 c2pv2'
    expression = '50.0*((c1pv2-0.187)^2+2*(c2pv2-0.0157)^2)'
  [../]
  # Elastic energy of the phase 2
  [./elastic_free_energy_pv2]
    type = ElasticEnergyMaterial
    base_name = phasepv2
    property_name = fe_pv2
    coupled_variables = ' '
  [../]
  # Total free energy of the phase 2
  [./Total_energy_pv2]
    type = DerivativeSumMaterial
    property_name = Fpv2
    sum_materials = 'fc_pv2 fe_pv2'
    coupled_variables = 'c1pv2 c2pv2'
  [../]

  [./f3]
    type = DerivativeParsedMaterial
    property_name = fc_pv3
    coupled_variables = 'c1pv3 c2pv3'
    expression = '50.0*((c1pv3-0.000727)^2+2*(c2pv3-0.196)^2)'
  [../]
    # Elastic energy of the phase 3
  [./elastic_free_energy_pv3]
    type = ElasticEnergyMaterial
    base_name = phasepv3
    property_name = fe_pv3
    coupled_variables = ' '
  [../]
    # Total free energy of the phase 3
  [./Total_energy_pv3]
    type = DerivativeSumMaterial
    property_name = Fpv3
    sum_materials = 'fc_pv3 fe_pv3'
    coupled_variables = 'c1pv3 c2pv3'
  [../]
   # --- δ chemical free energy (gated ON at ~3 h) ---
  # Lowered prefactor (35 < 50) makes δ thermodynamically more favorable than γ″,
  # matching that δ is the equilibrium phase. 'gate_3h' keeps it ~0 before 3 h.
  [./fc_pv4]
  type = DerivativeParsedMaterial
  property_name     = fc_pv4
  coupled_variables = 'c1pv4 c2pv4'
  # pick centers near your actual compositions (example below)
  expression        = '35*((c1pv4 - 0.020)^2 + 2*(c2pv4 - 0.060)^2)'
[../]

  # Elastic energy of the phase 4
  [./elastic_free_energy_pv4]
    type = ElasticEnergyMaterial
    base_name = phasepv4
    property_name = fe_pv4
    coupled_variables = ' '
  [../]
    # Total free energy of the phase 4
  [./Total_energy_pv4]
  type = DerivativeSumMaterial
  property_name  = Fpv4
  sum_materials  = 'fc_pv4 fe_pv4'
  coupled_variables = 'c1pv4 c2pv4'
[../]

  # Switching functions for each phase
  # hm(eta_pv1, eta_pv2, eta_m)
  [./hm]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_m
    all_etas = 'eta_pv1 eta_pv2 eta_pv3 eta_m'
    h_name = hm
  [../]
  # hpv1(eta_pv1, eta_pv2, eta_m)
  [./hpv1]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv1
    all_etas = 'eta_pv1 eta_pv2 eta_pv3 eta_m'
    h_name = hpv1
  [../]
  # hpv2(eta_pv1, eta_pv2, eta_m)
  [./hpv2]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv2
    all_etas = 'eta_pv1 eta_pv2 eta_pv3 eta_m'
    h_name = hpv2
  [../]
[./hpv3]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv3
    all_etas = 'eta_pv1 eta_pv2 eta_pv3 eta_m'
    h_name = hpv3
[../]
  # hpv4(eta_pv1, eta_pv2, eta_m)
  [./hpv4]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv4
    all_etas = 'eta_pv1 eta_pv2 eta_pv3 eta_m'
    h_name = hpv4
  [../]
  # Coefficients for diffusion equation
  [./Dhm]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D hm'
    expression = gb_scale_aux*(D*hm)
    property_name = Dhm
  [../]
  [./Dhpv1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D hpv1'
    expression = gb_scale_aux*(D*hpv1)
    property_name = Dhpv1
  [../]
  [./Dhpv2]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D hpv2'
    expression = gb_scale_aux*(D*hpv2)
    property_name = Dhpv2
  [../]
  [./Dhpv3]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D hpv3'
    expression = gb_scale_aux*(D*hpv3)
    property_name = Dhpv3
  [../]
  [./Dhpv4]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D hpv4'
    expression = gb_scale_aux*(D*hpv4)
    property_name = Dhpv4
  [../]

# Barrier functions for each phase
  [./gm]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta_m
    function_name = gm
  [../]
  [./gpv1]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_pv1
    function_name = gpv1
  [../]
  [./gpv2]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_pv2
    function_name = gpv2
  [../]
  [./gpv3]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_pv3
    function_name = gpv3
  [../]
  [./gpv4]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_pv4
    function_name = gpv4
  [../]

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'L    kappa kappa_pv4 D  misfit     W'
    prop_values = '0.3  0.01    0.008   1    1        0.01'
  [../]
  #[./gate_3h_prop]
   # type        = GenericFunctionMaterial
  #  prop_names  = 'gate_3h'
  #  prop_values = 'gate_3h'
  #[../]

  #Mechanical properties
  [./Stiffness_phasem_g0]
    type = ComputeElasticityTensor
    C_ijkl = '272.1 169 169 272.1 169 272.1 131 131 131' #Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of     
    base_name = phasem
    fill_method = symmetric9
    euler_angle_1 = 0   
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 0
  [../]
  [./Stiffness_phasem_g1]
    type = ComputeElasticityTensor
    C_ijkl = '272.1 169 169 272.1 169 272.1 131 131 131' #Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of     
    base_name = phasem
    fill_method = symmetric9
    euler_angle_1 = 45
    euler_angle_2  = 0
    euler_angle_3  = 0
    block = 1
  [../]
  [./Stiffness_phasepv1_g0]
  type = ComputeElasticityTensor
  base_name   = phasepv1
  C_ijkl      = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
  fill_method = symmetric9
  euler_angle_1 = 0
  euler_angle_2 = 0
  euler_angle_3 = 0
  block = 0
[../]
[./Stiffness_phasepv1_g1]
  type = ComputeElasticityTensor
  base_name   = phasepv1
  C_ijkl      = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
  fill_method = symmetric9
  euler_angle_1 = 45
  euler_angle_2 = 0
  euler_angle_3 = 0
  block = 1
[../]
[./Stiffness_phasepv2_g0]
  type = ComputeElasticityTensor
  C_ijkl = '243 154.8 154.8 243 154.8 243 132.3 132.3 132.3'
  base_name = phasepv2
  fill_method = symmetric9
  euler_angle_1 = 0
  euler_angle_2 = 0
  euler_angle_3 = 0
  block = 0
  [../]
  [./Stiffness_phasepv2_g1]
    type = ComputeElasticityTensor
    C_ijkl = '243 154.8 154.8 243 154.8 243 132.3 132.3 132.3'
    base_name = phasepv2
    fill_method = symmetric9
    euler_angle_1 = 45
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 1
  [../]
  [./Stiffness_phasepv3_g0]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
    base_name = phasepv3
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 0
  [../]
  [./Stiffness_phasepv3_g1]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
    base_name = phasepv3
    fill_method = symmetric9
    euler_angle_1 = 45
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 1
  [../]
  [./Stiffness_phasepv4_g0]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
    base_name = phasepv4
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
  [../]
  

  [./stress_phasepv1_g0]
    type = ComputeLinearElasticStress
    base_name = phasepv1 
    block = 0
  [../]
  [./stress_phasepv1_g1]
    type = ComputeLinearElasticStress
    base_name = phasepv1
    block = 1
  [../]
  [./stress_phasepv2_g0]
    type = ComputeLinearElasticStress
    base_name = phasepv2
    block = 0
  [../]
  [./stress_phasepv2_g1]
    type = ComputeLinearElasticStress
    base_name = phasepv2
    block = 1
  [../]
  [./stress_phasepv3_g0]
    type = ComputeLinearElasticStress
    base_name = phasepv3
    block = 0
  [../]
  [./stress_phasepv3_g1]
    type = ComputeLinearElasticStress
    base_name = phasepv3
    block = 1
  [../]
  [./stress_phasepv4_g0]
    type = ComputeLinearElasticStress
    base_name = phasepv4
  [../]

  [./stress_phasem_g0]
    type = ComputeLinearElasticStress
    base_name = phasem
    block = 0
  [../]
  [./stress_phasem_g1]
    type = ComputeLinearElasticStress
    base_name = phasem
    block = 1
  [../]

  [./strain_phasem_g0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasem
    block = 0
  [../]
  [./strain_phasem_g1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasem
    block = 1
  [../]
  [./strain_phasepv1_g0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv1
    eigenstrain_names = eigenstrainpv1
    block = 0
  [../]
  [./strain_phasepv1_g1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv1
    eigenstrain_names = eigenstrainpv1
    block = 1
  [../]
  [./strain_phasepv2_g0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv2
    eigenstrain_names = eigenstrainpv2
    block = 0
  [../]
  [./strain_phasepv2_g1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv2
    eigenstrain_names = eigenstrainpv2
    block = 1
  [../]
  [./strain_phasepv3_g0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv3
    eigenstrain_names = eigenstrainpv3
    block = 0
  [../]
  [./strain_phasepv3_g1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv3
    eigenstrain_names = eigenstrainpv3
    block = 1
  [../]
  [./strain_phasepv4_g0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv4
    eigenstrain_names = eigenstrainpv4
  [../]
 



  [./eigen_strainpv1_g0]
    type = ComputeRotatedEigenstrain
    base_name = phasepv1
    eigen_base = '0.028 0.0067 0 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv1
    block = 0
  [../]
  [./eigen_strainpv1_g1]
    type = ComputeRotatedEigenstrain
    base_name = phasepv1
    eigen_base = '0.028 0.0067 0 0 0 0'
    Euler_angles = '45 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv1
    block = 1
  [../]

  [./eigen_strainpv2_g0]
    type = ComputeRotatedEigenstrain
    base_name = phasepv2
    eigen_base = '-0.003 -0.003 0 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv2
    block = 0
  [../]
   [./eigen_strainpv2_g1]
    type = ComputeRotatedEigenstrain
    base_name = phasepv2
    eigen_base = '-0.003 -0.003 0 0 0 0'
    Euler_angles = '45 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv2
    block = 1
  [../]
  [./eigen_strainpv3_g0]
    type = ComputeRotatedEigenstrain
    base_name = phasepv3
    eigen_base = '0.0067 0.028 0 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv3
    block = 0
  [../]
  [./eigen_strainpv3_g1]
    type = ComputeRotatedEigenstrain
    base_name = phasepv3
    eigen_base = '0.0067 0.028 0 0 0 0'
    Euler_angles = '45 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv3
    block = 1
  [../]
  [./eigen_strainpv4_g0]
    type = ComputeRotatedEigenstrain
    base_name = phasepv4
    eigen_base = '0.005320 0.01332 0.02378 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = 1
    eigenstrain_name = eigenstrainpv4
  [../]
  


  # Generate the global stress from the phase stresses
  [./global_stress]
    type = MultiPhaseStressMaterial
    phase_base = 'phasepv1 phasepv2 phasepv3 phasepv4 phasem'
    h          = 'hpv1     hpv2   hpv3   hpv4   hm'
  [../]

  [./global_strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
  
  # Kernels for Allen-Cahn equation for eta_pv1
  [./deta_pv1_dt]
    type = TimeDerivative
    variable = eta_pv1
  [../]
  [./ACBulkFpv1]
    type = KKSMultiACBulkF
    variable  = eta_pv1
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    gi_name   = gpv1
    eta_i     = eta_pv1
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1pv4 c1m c2pv1 c2pv2 c2pv3 c2pv4 c2m eta_pv2 eta_pv3 eta_pv4 eta_m'
  [../]
  [./ACBulkCpv1_c1]
    type = KKSMultiACBulkC
    variable  = eta_pv1
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1pv4 c1m'
    eta_i     = eta_pv1
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2pv4 c2m eta_pv2 eta_pv3 eta_pv4 eta_m'
  [../]
  [./ACBulkCpv1_c2]
    type = KKSMultiACBulkC
    variable  = eta_pv1
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2pv4 c2m'
    eta_i     = eta_pv1
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1pv4 c1m  eta_pv2 eta_pv3 eta_pv4 eta_m'
  [../]
  [./ACInterfacepv1]
    type = ACInterface
    variable = eta_pv1
    kappa_name = kappa
  [../]

  # Kernels for Allen-Cahn equation for eta_pv2
  [./deta_pv2_dt]
    type = TimeDerivative
    variable = eta_pv2
  [../]
  [./ACBulkFpv2]
    type = KKSMultiACBulkF
    variable  = eta_pv2
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    gi_name   = gpv2
    eta_i     = eta_pv2
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1pv4 c1m c2pv1 c2pv2 c2pv3 c2pv4 c2m eta_pv1 eta_pv3 eta_pv4 eta_m'
  [../]
  [./ACBulkCpv2_c1]
    type = KKSMultiACBulkC
    variable  = eta_pv2
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1pv4 c1m'
    eta_i     = eta_pv2
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2pv4 c2m eta_pv1 eta_pv3 eta_pv4 eta_m'
  [../]
  [./ACBulkCpv2_c2]
    type = KKSMultiACBulkC
    variable  = eta_pv2
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2pv4 c2m'
    eta_i     = eta_pv2
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1pv4 c1m eta_pv1 eta_pv3 eta_pv4 eta_m'
  [../]
  [./ACInterfacepv2]
    type = ACInterface
    variable = eta_pv2
    kappa_name = kappa
  [../]

  # Kernels for Allen-Cahn equation for eta_pv3
  [./deta_pv3_dt]
    type = TimeDerivative
    variable = eta_pv3
  [../]
  [./ACBulkFpv3]
    type = KKSMultiACBulkF
    variable  = eta_pv3
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    gi_name   = gpv3
    eta_i     = eta_pv3
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1pv4 c1m c2pv1 c2pv2 c2pv3 c2pv4 c2m eta_pv1 eta_pv2 eta_pv4 eta_m'
  [../]
  [./ACBulkCpv3_c1]
    type = KKSMultiACBulkC
    variable  = eta_pv3
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1pv4 c1m'
    eta_i     = eta_pv3
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2pv4 c2m eta_pv1 eta_pv2 eta_pv4 eta_m'
  [../]
  [./ACBulkCpv3_c2]
    type = KKSMultiACBulkC
    variable  = eta_pv3
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2pv4 c2m'
    eta_i     = eta_pv3
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1pv4 c1m eta_pv1 eta_pv2 eta_pv4 eta_m'
  [../]
  [./ACInterfacepv3]
    type = ACInterface
    variable = eta_pv3
    kappa_name = kappa
  [../]
 
# Kernels for Allen-Cahn equation for eta_pv4
  [./deta_pv4_dt]
    type = TimeDerivative
    variable = eta_pv4
  [../]
  [./ACBulkFpv4]
    type = KKSMultiACBulkF
    variable  = eta_pv4
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    gi_name   = gpv4
    eta_i     = eta_pv4
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1pv4 c1m c2pv1 c2pv2 c2pv3 c2pv4 c2m eta_pv1 eta_pv2 eta_pv3 eta_m'
  [../]
  [./ACBulkCpv4_c1]
    type = KKSMultiACBulkC
    variable  = eta_pv4
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1pv4 c1m'
    eta_i     = eta_pv4
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2pv4 c2m eta_pv1 eta_pv2 eta_pv3 eta_m'
  [../]
  [./ACBulkCpv4_c2]
    type = KKSMultiACBulkC
    variable  = eta_pv4
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2pv4 c2m'
    eta_i     = eta_pv4
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1pv4 c1m eta_pv1 eta_pv2 eta_pv3 eta_m'
  [../]
  [./ACInterfacepv4]
    type = ACInterface
    variable = eta_pv4
    kappa_name = kappa_pv4
  [../]

# Kernels for constraint equation |eta_pv1| + |eta_pv2| + eta_m = 1
  # eta3 is the nonlinear variable for the constraint equation
  [./eta_mreaction]
    type = MatReaction
    variable = eta_m
    reaction_rate = 1
 [../]
  [./eta_pv1reaction]
    type = MatReaction_abscouple
    variable = eta_m
    v = eta_pv1
    reaction_rate = 1
  [../]
  [./eta_pv2reaction]
    type = MatReaction_abscouple
    variable = eta_m
    v = eta_pv2
    reaction_rate = 1
  [../]
  [./eta_pv3reaction]
    type = MatReaction_abscouple
    variable = eta_m
    v = eta_pv3
    reaction_rate = 1
  [../]
  [./eta_pv4reaction]
    type = MatReaction_abscouple
    variable = eta_m
    v = eta_pv4
    reaction_rate = 1
  [../]
  [./on]
    type = BodyForce
    variable = eta_m
    value = -1.0
  [../]

  #Kernels for diffusion equation of c1
  [./diff_time_c1]
    type = TimeDerivative
    variable = c1
  [../]
  [./diff_c1m]
    type = MatDiffusion
    variable = c1
    diffusivity = Dhm
    v = c1m
  [../]
  [./diff_c1pv1]
    type = MatDiffusion
    variable = c1
    diffusivity = Dhpv1
    v = c1pv1
  [../]
  [./diff_c1pv2]
    type = MatDiffusion
    variable = c1
    diffusivity = Dhpv2
    v = c1pv2
  [../]
  [./diff_c1pv3]
    type = MatDiffusion
    variable = c1
    diffusivity = Dhpv3
    v = c1pv3
  [../]
  [./diff_c1pv4]
    type = MatDiffusion
    variable = c1
    diffusivity = Dhpv4
    v = c1pv4
  [../]


  #Kernels for diffusion equation of c2
  [./diff_time_c2]
    type = TimeDerivative
    variable = c2
  [../]
  [./diff_c2m]
    type = MatDiffusion
    variable = c2
    diffusivity = Dhm
    v = c2m
  [../]
  [./diff_c2pv1]
    type = MatDiffusion
    variable = c2
    diffusivity = Dhpv1
    v = c2pv1
  [../]
  [./diff_c2pv2]
    type = MatDiffusion
    variable = c2
    diffusivity = Dhpv2
    v = c2pv2
  [../]
  [./diff_c2pv3]
    type = MatDiffusion
    variable = c2
    diffusivity = Dhpv3
    v = c2pv3
  [../]   
  [./diff_c2pv4]
    type = MatDiffusion
    variable = c2
    diffusivity = Dhpv4
    v = c2pv4
  [../]

  # Phase concentration constraints
   [./chempot1m_pv1]
    type = KKSPhaseChemicalPotential
    variable = c1m
    cb       = c1pv1
    fa_name  = Fm
    fb_name  = Fpv1
    args_a   = c2m
    args_b   = c2pv1
  [../]
 [./chempot1m_pv2]
    type = KKSPhaseChemicalPotential
    variable = c1pv1
    cb       = c1pv2
    fa_name  = Fpv1
    fb_name  = Fpv2
    args_a   = c2pv1
    args_b   = c2pv2
  [../]
  [./chempot1m_pv3]
    type = KKSPhaseChemicalPotential
    variable = c1pv2
    cb       = c1pv3
    fa_name  = Fpv2
    fb_name  = Fpv3
    args_a   = c2pv2
    args_b   = c2pv3
  [../]
  [./chempot1m_pv4]
    type = KKSPhaseChemicalPotential
    variable = c1pv3
    cb       = c1pv4
    fa_name  = Fpv3
    fb_name  = Fpv4
    args_a   = c2pv3
    args_b   = c2pv4
  [../]
  [./chempot2m_pv1]
    type = KKSPhaseChemicalPotential
    variable = c2m
    cb       = c2pv1
    fa_name  = Fm
    fb_name  = Fpv1
    args_a   = c1m
    args_b   = c1pv1
  [../]
  [./chempot2m_pv2]
    type = KKSPhaseChemicalPotential
    variable = c2pv1
    cb       = c2pv2
    fa_name  = Fpv1
    fb_name  = Fpv2
    args_a   = c1pv1
    args_b   = c1pv2
  [../]
  [./chempot2m_pv3]
    type = KKSPhaseChemicalPotential
    variable = c2pv2
    cb       = c2pv3
    fa_name  = Fpv2
    fb_name  = Fpv3
    args_a   = c1pv2
    args_b   = c1pv3
  [../]
  [./chempot2m_pv4]
    type = KKSPhaseChemicalPotential
    variable = c2pv3
    cb       = c2pv4
    fa_name  = Fpv3
    fb_name  = Fpv4
    args_a   = c1pv3
    args_b   = c1pv4
  [../]
    
  [./phaseconcentration_c1pv3]
    type = KKSMultiPhaseConcentration
    variable = c1pv3
    cj = 'c1m c1pv1 c1pv2 c1pv3 c1pv4'
    hj_names = 'hm hpv1 hpv2 hpv3 hpv4'
    etas = 'eta_m eta_pv1 eta_pv2 eta_pv3 eta_pv4'
    c = c1
  [../]

  [./phaseconcentration_c1pv4]
    type = KKSMultiPhaseConcentration
    variable = c1pv4
    cj = 'c1m c1pv1 c1pv2 c1pv3 c1pv4'
    hj_names = 'hm hpv1 hpv2 hpv3 hpv4'
    etas = 'eta_m eta_pv1 eta_pv2 eta_pv3 eta_pv4'
    c = c1
  [../]

  [./phaseconcentration_c2pv3]
    type = KKSMultiPhaseConcentration
    variable = c2pv3
    cj = 'c2m c2pv1 c2pv2 c2pv3 c2pv4'
    hj_names = 'hm hpv1 hpv2 hpv3 hpv4'
    etas = 'eta_m eta_pv1 eta_pv2 eta_pv3 eta_pv4'
    c = c2
  [../]
    
  [./phaseconcentration_c2pv4]
    type = KKSMultiPhaseConcentration
    variable = c2pv4
    cj = 'c2m c2pv1 c2pv2 c2pv3 c2pv4'
    hj_names = 'hm hpv1 hpv2 hpv3 hpv4'
    etas = 'eta_m eta_pv1 eta_pv2 eta_pv3 eta_pv4'
    c = c2
  [../]
  

[]

[AuxKernels]
[./c2_gb_enrich_fn]
    type      = FunctionAux
    variable  = gb_mask
    function  = c2_gb_enrich_fn
  [../]

  # 2) Weighted c2 for integration/averaging over GB band
  [./c2_times_mask]
    type               = ParsedAux
    variable           = c2_gb_w
    coupled_variables  = 'c2 gb_mask'
    expression         = 'c2 * gb_mask'
  [../]

  # 3) Weighted ∆c2 relative to a baseline (set baseline to your far-field Nb, e.g. 0.038)
  [./dc2_times_mask]
    type               = ParsedAux
    variable           = dc2_gb_w
    coupled_variables  = 'c2 gb_mask'
    expression         = '(c2 - 0.038) * gb_mask'
  [../]
  [./gb_scale_eval]
    type     = FunctionAux
    variable = gb_scale_aux
    function = gb_scale_fn
  [../]
  [./set_one] 
   type = ParsedAux 
   variable = one 
   expression = '1' 
   execute_on = timestep_end 
  [../]
    [./copy_g]
    type = FunctionAux
    variable = g_aux
    function = gb_scale_fn
    execute_on = timestep_end
  [../]
# weighted hpv by g (use your actual hpv variable names)
  [./mk_h_g_pv1] 
  type=ParsedAux 
  variable=h_g_pv1 
  coupled_variables = 'g_aux h_pv1_aux'
  expression='g_aux * h_pv1_aux' 
  execute_on=timestep_end 
  [../]
  [./mk_h_g_pv2] 
  type=ParsedAux 
  variable=h_g_pv2 
  coupled_variables = 'g_aux h_pv2_aux'
  expression='g_aux * h_pv2_aux' 
  execute_on=timestep_end 
  [../]
  [./mk_h_g_pv3] 
  type=ParsedAux 
  variable=h_g_pv3 
  coupled_variables = 'g_aux h_pv3_aux'
  expression='g_aux * h_pv3_aux' 
  execute_on=timestep_end 
  [../]
  [./mk_h_g_pv4]
  type=ParsedAux
  variable=h_g_pv4
  coupled_variables = 'g_aux h_pv4_aux'
  expression='g_aux * h_pv4_aux'
  execute_on=timestep_end
  [../]
  [temperature]
    type = FunctionAux
    variable = temperature
    function = '1073'
    execute_on = timestep_begin
  []
  [./Energy_total]
    type = KKSMultiFreeEnergy
    Fj_names = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names = 'hpv1 hpv2 hpv3 hpv4 hm'
    gj_names = 'gpv1 gpv2 gpv3 gpv4 gm'
    variable = Energy
    w = 1
    interfacial_vars =  'eta_pv1  eta_pv2  eta_pv3  eta_pv4  eta_m'
    kappa_names =       'kappa kappa kappa kappa kappa'
  [../]
   [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    scalar_type = VonMisesStress
    execute_on = timestep_end
  [../]
   [./copy_h_pv1] 
  type = MaterialRealAux
  variable = h_pv1_aux 
  property = hpv1 
  execute_on = timestep_end 
  [../]
  [./copy_h_pv2] 
  type = MaterialRealAux 
  variable = h_pv2_aux 
  property = hpv2 
  execute_on = timestep_end 
  [../]
  [./copy_h_pv3] 
  type = MaterialRealAux 
  variable = h_pv3_aux 
  property = hpv3 
  execute_on = timestep_end 
  [../]
  [./copy_h_pv4]
  type = MaterialRealAux
  variable = h_pv4_aux
  property = hpv4
  execute_on = timestep_end
  [../]
  [./copy_h_m]
  type = MaterialRealAux
  variable = h_m_aux
  property = hm
  execute_on = timestep_end
  [../]
   # Products: h * σ_xx
   [./vonmises_times_h_pv1]
     type = ParsedAux
     variable = vonmises_h_pv1
     coupled_variables = 'vonmises h_pv1_aux'
     expression = 'vonmises * h_pv1_aux'
     execute_on = timestep_end
   [../]
  [./vonmises_times_h_pv2]
      type = ParsedAux
      variable = vonmises_h_pv2
      coupled_variables = 'vonmises h_pv2_aux'
      expression = 'vonmises * h_pv2_aux'
      execute_on = timestep_end
    [../]
    [./vonmises_times_h_pv3]
      type = ParsedAux
      variable = vonmises_h_pv3
      coupled_variables = 'vonmises h_pv3_aux'
      expression = 'vonmises * h_pv3_aux'
      execute_on = timestep_end
    [../]
    [./vonmises_times_h_pv4]
      type = ParsedAux
      variable = vonmises_h_pv4
      coupled_variables = 'vonmises h_pv4_aux'
      expression = 'vonmises * h_pv4_aux'
      execute_on = timestep_end
    [../]
    [./vonmises_times_h_m]
      type = ParsedAux
      variable = vonmises_h_m
      coupled_variables = 'vonmises h_m_aux'
      expression = 'vonmises * h_m_aux'
      execute_on = timestep_end
    [../]
  
  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
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
  [./mydebug]
    type = FDP
    full = true
  [../]
[]

[Postprocessors]
   [./dofs]
     type = NumDOFs
   [../]
   [./h1_error]
     type = ElementH1Error
     variable = eta_pv1
     function = bc_func
   [../]
  [./Energy_total]
    type = ElementAverageValue
    variable = Energy
  []
  [stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  []
   [./vonmises]
    type = ElementAverageValue
    variable = vonmises
  [../]
  [./num_vonmises_pv1]
  type = ElementIntegralVariablePostprocessor
  variable = vonmises_h_pv1
  [../]
  [./num_vonmises_pv2]
  type = ElementIntegralVariablePostprocessor
  variable = vonmises_h_pv2
  [../]
  [./num_vonmises_pv3]
  type = ElementIntegralVariablePostprocessor
  variable = vonmises_h_pv3
  [../]
  [./num_vonmises_pv4]
  type = ElementIntegralVariablePostprocessor
  variable = vonmises_h_pv4
  [../]
  [./num_vonmises_m]  
   type = ElementIntegralVariablePostprocessor
   variable = vonmises_h_m   
   [../]
  

  # Denominators: ∫ h dV  (phase volumes)
  [./den_pv1] 
  type = ElementIntegralVariablePostprocessor 
  variable = h_pv1_aux 
  [../]
  [./den_pv2] 
  type = ElementIntegralVariablePostprocessor 
  variable = h_pv2_aux 
  [../]
  [./den_pv3] 
  type = ElementIntegralVariablePostprocessor 
  variable = h_pv3_aux 
  [../]
  [./den_pv4] 
  type = ElementIntegralVariablePostprocessor 
  variable = h_pv4_aux 
  [../]
  [./den_m]   
   type = ElementIntegralVariablePostprocessor 
   variable = h_m_aux   
   [../]
  [./eta_pv1]
    type = ElementIntegralVariablePostprocessor
    variable = eta_pv1
    use_absolute_value = true
  [../]
   [./eta_pv2]
    type = ElementIntegralVariablePostprocessor
    variable = eta_pv2
    use_absolute_value = true
  [../]
  [./eta_pv3]
    type = ElementIntegralVariablePostprocessor
    variable = eta_pv3
    use_absolute_value = true
  [../]
  [./eta_pv4]
    type = ElementIntegralVariablePostprocessor
    variable = eta_pv4
    use_absolute_value = true
  [../]
   [./A_total]   
  type=ElementIntegralVariablePostprocessor 
  variable=one       
  [../]
  [./A_g]       
  type=ElementIntegralVariablePostprocessor 
  variable=g_aux     
  [../]
    # Average hpv in whole domain
    [./af_pv1] 
      type = ElementAverageMaterialProperty 
      mat_prop = hpv1 
    [../]
    [./af_pv2] 
      type = ElementAverageMaterialProperty 
      mat_prop = hpv2 
    [../]
    [./af_pv3] 
      type = ElementAverageMaterialProperty 
      mat_prop = hpv3 
    [../]
  [./af_pv4]
    type = ElementAverageMaterialProperty
    mat_prop = hpv4
  [../]
  [./den_g_pv1] 
  type=ElementIntegralVariablePostprocessor 
  variable=h_g_pv1   
  [../]
  [./den_g_pv2] 
  type=ElementIntegralVariablePostprocessor 
  variable=h_g_pv2   
  [../]
  [./den_g_pv3] 
  type=ElementIntegralVariablePostprocessor 
  variable=h_g_pv3   
  [../]
  [./den_g_pv4]
  type=ElementIntegralVariablePostprocessor
  variable=h_g_pv4
  [../]
[./int_c2_func]
    type      = FunctionElementIntegral
    function  = c2_gb_enrich_fn
  [../]
   # fractions in GB band and in whole domain
  [./frac_g_pv1] 
  type=ParsedPostprocessor 
  pp_names='den_g_pv1 A_g'      
  expression='den_g_pv1/max(A_g,1e-16)'     
  [../]
  [./frac_g_pv2] 
  type=ParsedPostprocessor 
  pp_names='den_g_pv2 A_g'     
   expression='den_g_pv2/max(A_g,1e-16)'     
   [../]
  [./frac_g_pv3] 
  type=ParsedPostprocessor 
  pp_names='den_g_pv3 A_g'      
  expression='den_g_pv3/max(A_g,1e-16)'     
  [../]
  [./frac_pv1]   
  type=ParsedPostprocessor 
  pp_names='den_pv1 A_total'     
  expression='den_pv1/max(A_total,1e-16)'   
  [../]
  [./frac_pv2]   
  type=ParsedPostprocessor 
  pp_names='den_pv2 A_total'     
  expression='den_pv2/max(A_total,1e-16)'   
  [../]
  [./frac_pv3]   
  type=ParsedPostprocessor 
  pp_names='den_pv3 A_total'     
  expression='den_pv3/max(A_total,1e-16)'  
   [../]
  [./frac_pv4]   
  type=ParsedPostprocessor 
  pp_names='den_pv4 A_total'     
  expression='den_pv4/max(A_total,1e-16)'  
   [../]
[]

[Outputs]
  exodus = true
  [./table]
    type = CSV
    execute_on = timestep_end
  [../]
  [./checkpoint]
     type = Checkpoint
     num_files = 10
     time_step_interval = 10
  [../]
[]
