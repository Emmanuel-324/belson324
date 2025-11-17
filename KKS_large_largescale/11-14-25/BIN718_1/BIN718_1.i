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
  [./eta_gp_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_gp
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_gp_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_gp
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_gpp1_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_gpp1
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_gpp1_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_gpp1
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_gpp2_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_gpp2
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_gpp2_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_gpp2
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_d_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_d
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_d_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_d
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
# phase concentration c1 in gp
  [./c1gp]
    order = FIRST
    family = LAGRANGE
  [../]
  # phase concentration c1 in gpp1
  [./c1gpp1]
    order = FIRST
    family = LAGRANGE
  [../]
 # phase concentration c1 in gpp2
 [./c1gpp2]
  order = FIRST
  family = LAGRANGE
[../]
  # phase concentration c1 in d
  [./c1d]
    order = FIRST
    family = LAGRANGE
  [../]

# phase concentration c2 in matrix
  [./c2m]
    order = FIRST
    family = LAGRANGE
  [../]
  # phase concentration c2 in gp
  [./c2gp]
    order = FIRST
    family = LAGRANGE
  [../]
  # phase concentration c2 in gpp1
  [./c2gpp1]
    order = FIRST
    family = LAGRANGE
  [../]
 # phase concentration c2 in gpp2
 [./c2gpp2]
  order = FIRST
  family = LAGRANGE
[../]
  # phase concentration c2 in d
  [./c2d]
    order = FIRST
    family = LAGRANGE
  [../]

  # order parameter m
  [./eta_m]
    order = FIRST
    family = LAGRANGE
  [../]
  # order parameter pv1
  [./eta_gp]
    order = FIRST
    family = LAGRANGE
  [../]
  # order parameter pv2
  [./eta_gpp1]
    order = FIRST
    family = LAGRANGE
  [../]
# order parameter pv3
[./eta_gpp2]
  order = FIRST
  family = LAGRANGE
[../]
  # order parameter pv4
  [./eta_d]
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
    value = '0.05*exp(-((x - 350)/6.0)^2)'            # ~5 nm half-width, amplitude 0.02
  [../] 
# ---- Nb (c2) enrichment at the grain boundary (x = 350 sim units) ----
  # Baseline values from Wang: x_Al ≈ 0.024, x_Nb ≈ 0.038 (adjust if your alloy differs).
  # dC_nb is the peak enrichment (keep modest to preserve mass balance).
 
[]

[ICs]
   [./eta_gp]
    variable = eta_gp
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 324
  [../]
  [./eta_gpp1]
    variable = eta_gpp1
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 230	
  [../]
  [./eta_gpp2]
    variable = eta_gpp2
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 307	
  [../]
  # thin delta plate centered on the GB (spans both sides a little)
[./eta_d_gb_plate_ic]
    type     = FunctionIC
    variable = eta_d
    function = '0.8 * 0.5*(tanh((x - 340.0)/1.0) - tanh((x - 360.0)/1.0))'
  [../]
  [./c1]
    variable = c1
    type = RandomIC
    min = 0.010	
    max = 0.03
    seed = 403	
  [../]
  [./c2]
    variable = c2
    type = RandomIC
    min = 0.032	
    max = 0.044
    seed = 89	
  [../]
[]


[Materials]
 # ------------------------------------------------------------------
  # 1.  Chemical free energy of the matrix (γ)  
  # ------------------------------------------------------------------
  [./fc_gamma]
    type = DerivativeParsedMaterial
    property_name = fc_gamma
    coupled_variables = 'c1m c2m'               # c1m = Al, c2m = Nb  (γ-phase)
    constant_names = 'Vm_inv C_Al C_Nb xeq_g_Al xeq_g_Nb'
    constant_expressions = '1e5 4.0e5 9.5e5 0.0161 0.00723'
    expression = 'Vm_inv*(C_Al*(c1m-xeq_g_Al)^2 + C_Nb*(c2m-xeq_g_Nb)^2)'
    derivative_order = 2
  [../]
  # Elastic energy of the phase 0
  [./elastic_free_energy_m]
    type = ElasticEnergyMaterial
    base_name = phasem
    property_name = fe_gamma
    coupled_variables = ' '
  [../]
  # Total free energy of the phase 0
  [./Total_energy_m]
    type = DerivativeSumMaterial
    property_name = F_gamma
    sum_materials = 'fc_gamma fe_gamma'
    coupled_variables = 'c1m c2m'
  [../]
# ------------------------------------------------------------------
  # 2.  Chemical free energy of γ' (L1₂)  
  # ------------------------------------------------------------------
  [./fc_gp]
    type = DerivativeParsedMaterial
    property_name = fc_gp
    coupled_variables = 'c1gp c2gp'             # γ'-phase concentrations
    constant_names = 'Vm_inv C_Al C_Nb xeq_gp_Al xeq_gp_Nb'
    constant_expressions = '1e5 4.0e5 9.5e5 0.187 0.0157'
    expression = 'Vm_inv*(C_Al*(c1gp-xeq_gp_Al)^2 + C_Nb*(c2gp-xeq_gp_Nb)^2)'
    derivative_order = 2
  [../]
  # Elastic energy of the phase 1
  [./elastic_free_energy_gp]
    type = ElasticEnergyMaterial
    base_name = phasegp
    property_name = fe_gp
    coupled_variables = ' '
  [../]
  # Total free energy of the phase 1
  [./Total_energy_gp]
    type = DerivativeSumMaterial
    property_name = F_gp
    sum_materials = 'fc_gp fe_gp'
    coupled_variables = 'c1gp c2gp'
  [../]
# ------------------------------------------------------------------
  # 3.  γ'' variant 1 (metastable) = δ + Δf   
  # ------------------------------------------------------------------
 [./fc_gpp1]
    type = DerivativeParsedMaterial
    property_name = fc_gpp1
    coupled_variables = 'c1gpp1 c2gpp1'      # Al and Nb in γ'' phase
    constant_names = 'Vm_inv C_Al C_Nb xeq_d_Al xeq_d_Nb Delta_f'
    constant_expressions = '1e5 4.0e5 9.5e5 7.27e-4 0.196 1548.5'
    expression = ' Vm_inv * (C_Al*(c1gpp1 - xeq_d_Al)^2 + C_Nb*(c2gpp1 - xeq_d_Nb)^2) + Vm_inv * Delta_f'
    derivative_order = 2
  [../] 
  # Elastic energy of the phase 2
  [./elastic_free_energy_gpp1]
    type = ElasticEnergyMaterial
    base_name = phasegpp1
    property_name = fe_gpp1
    coupled_variables = ' '
  [../]
  # Total free energy of the phase 2
  [./Total_energy_gpp1]
    type = DerivativeSumMaterial
    property_name = F_gpp1
    sum_materials = 'fc_gpp1 fe_gpp1'
    coupled_variables = 'c1gpp1 c2gpp1'
  [../]
 
# ------------------------------------------------------------------
  # 4.  γ'' variant 2 (metastable) = δ + Δf   
  # ------------------------------------------------------------------
 [./fc_gpp2]
    type = DerivativeParsedMaterial
    property_name = fc_gpp2
    coupled_variables = 'c1gpp2 c2gpp2'      # Al and Nb in γ'' phase
    constant_names = 'Vm_inv C_Al C_Nb xeq_d_Al xeq_d_Nb Delta_f'
    constant_expressions = '1e5 4.0e5 9.5e5 7.27e-4 0.196 1548.5'
    expression = ' Vm_inv * (C_Al*(c1gpp2 - xeq_d_Al)^2 + C_Nb*(c2gpp2 - xeq_d_Nb)^2) + Vm_inv * Delta_f'
    derivative_order = 2
  [../] 
  # Elastic energy of the phase 2
  [./elastic_free_energy_gpp2]
    type = ElasticEnergyMaterial
    base_name = phasegpp2
    property_name = fe_gpp2
    coupled_variables = ' '
  [../]
  # Total free energy of the phase 2
  [./Total_energy_gpp2]
    type = DerivativeSumMaterial
    property_name = F_gpp2
    sum_materials = 'fc_gpp2 fe_gpp2'
    coupled_variables = 'c1gpp1 c2gpp1'
  [../]
# ------------------------------------------------------------------
  # 5.  Chemical free energy of δ (D0ₐ)  
  # ------------------------------------------------------------------  
[./fc_d]
    type = DerivativeParsedMaterial
    property_name = fc_d
    coupled_variables = 'c1d c2d'               # δ-phase concentrations
    constant_names = 'Vm_inv C_Al C_Nb xeq_d_Al xeq_d_Nb'
    constant_expressions = '1e5 4.0e5 9.5e5 7.27e-4 0.196'
    expression = 'Vm_inv*(C_Al*(c1d-xeq_d_Al)^2 + C_Nb*(c2d-xeq_d_Nb)^2)'
    derivative_order = 2
  [../]



  # Elastic energy of the phase 4
  [./elastic_free_energy_pv4]
    type = ElasticEnergyMaterial
    base_name = phased
    property_name = fe_d
    coupled_variables = ' '
  [../]
    # Total free energy of the phase 4
  [./Total_energy_pv4]
  type = DerivativeSumMaterial
  property_name  = F_d
  sum_materials  = 'fc_d fe_d'
  coupled_variables = 'c1d c2d'
[../]

  # Switching functions for each phase
  # hm(eta_gp, eta_gpp1, eta_m)
  [./hm]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_m
    all_etas = 'eta_gp eta_gpp1 eta_gpp2 eta_d eta_m'
    h_name = hm
  [../]
  # hgp(eta_gp, eta_gpp1, eta_m)
  [./hgp]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_gp
    all_etas = 'eta_gp eta_gpp1 eta_gpp2 eta_d eta_m'
    h_name = hgp
  [../]
  # hgpp1(eta_gp, eta_gpp1, eta_m)
  [./hgpp1]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_gpp1
    all_etas = 'eta_gp eta_gpp1 eta_gpp2 eta_d eta_m'
    h_name = hgpp1
  [../]
[./hgpp2]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_gpp2
    all_etas = 'eta_gp eta_gpp1 eta_gpp2 eta_d eta_m'
    h_name = hgpp2
[../]
  # hd(eta_gp, eta_gpp1, eta_m)
  [./hd]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_d
    all_etas = 'eta_gp eta_gpp1 eta_gpp2 eta_d eta_m'
    h_name = hd
  [../]
  # Coefficients for diffusion equation
  [./Dhm]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D hm'
    expression = gb_scale_aux*(D*hm)
    property_name = Dhm
  [../]
  [./Dhgp]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D hgp'
    expression = gb_scale_aux*(D*hgp)
    property_name = Dhgp
  [../]
  [./Dhgpp1]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D hgpp1'
    expression = gb_scale_aux*(D*hgpp1)
    property_name = Dhgpp1
  [../]
  [./Dhgpp2]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D hgpp2'
    expression = gb_scale_aux*(D*hgpp2)
    property_name = Dhgpp2
  [../]
  [./Dhd]
    type = DerivativeParsedMaterial
    coupled_variables = 'gb_scale_aux'
    material_property_names = 'D hd'
    expression = gb_scale_aux*(D*hd)
    property_name = Dhd
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
    eta = eta_gp
    function_name = gpv1
  [../]
  [./gpv2]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_gpp1
    function_name = gpv2
  [../]
  [./gpv3]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_gpp2
    function_name = gpv3
  [../]
  [./gpv4]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_d
    function_name = gpv4
  [../]

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'L    kappa  D   misfit    W'
    prop_values = '0.3   0.01  1     1      0.01'
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
  [./Stiffness_phasegp_g0]
  type = ComputeElasticityTensor
  C_ijkl = '243 154.8 154.8 243 154.8 243 132.3 132.3 132.3'
  base_name = phasegp
  fill_method = symmetric9
  euler_angle_1 = 0
  euler_angle_2 = 0
  euler_angle_3 = 0
  block = 0
  [../]
  [./Stiffness_phasegp_g1]
    type = ComputeElasticityTensor
    C_ijkl = '243 154.8 154.8 243 154.8 243 132.3 132.3 132.3'
    base_name = phasegp
    fill_method = symmetric9
    euler_angle_1 = 45
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 1
  [../]
  [./Stiffness_phasegpp1_g0]
  type = ComputeElasticityTensor
  base_name   = phasegpp1
  C_ijkl      = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
  fill_method = symmetric9
  euler_angle_1 = 0
  euler_angle_2 = 0
  euler_angle_3 = 0
  block = 0
[../]
[./Stiffness_phasegpp1_g1]
  type = ComputeElasticityTensor
  base_name   = phasegpp1
  C_ijkl      = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
  fill_method = symmetric9
  euler_angle_1 = 45
  euler_angle_2 = 0
  euler_angle_3 = 0
  block = 1
[../]

  [./Stiffness_phasegpp2_g0]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
    base_name = phasegpp2
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 0
  [../]
  [./Stiffness_phasegpp2_g1]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
    base_name = phasegpp2
    fill_method = symmetric9
    euler_angle_1 = 45
    euler_angle_2 = 0
    euler_angle_3 = 0
    block = 1
  [../]
   # per-block stiffness for pv4
[./Stiffness_phased_g0]
  type = ComputeElasticityTensor
  C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
  base_name = phased
  fill_method = symmetric9
  euler_angle_1 = 5.176036589
  euler_angle_2 = 1.150261992
  euler_angle_3 = 2.456873451
  block = 0
[../]

[./Stiffness_phased_g1]
  type = ComputeElasticityTensor
  C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
  base_name = phased
  fill_method = symmetric9
  euler_angle_1 = 5.176036589
  euler_angle_2 = 1.150261992
  euler_angle_3 = 2.456873451
  block = 1
[../]
  
[./stress_phasegp_g0]
    type = ComputeLinearElasticStress
    base_name = phasegp
    block = 0
  [../]
  [./stress_phasegp_g1]
    type = ComputeLinearElasticStress
    base_name = phasegp
    block = 1
  [../]
  [./stress_phasegpp1_g0]
    type = ComputeLinearElasticStress
    base_name = phasegpp1 
    block = 0
  [../]
  [./stress_phasegpp1_g1]
    type = ComputeLinearElasticStress
    base_name = phasegpp1
    block = 1
  [../]
  
  [./stress_phasegpp2_g0]
    type = ComputeLinearElasticStress
    base_name = phasegpp2
    block = 0
  [../]
  [./stress_phasegpp2_g1]
    type = ComputeLinearElasticStress
    base_name = phasegpp2
    block = 1
  [../]
   [./stress_phased_g0]
    type = ComputeLinearElasticStress
    base_name = phased
    block = 0
  [../]
  [./stress_phased_g1]
    type = ComputeLinearElasticStress
    base_name = phased
    block = 1
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
    [./strain_phasegp_g0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasegp
    eigenstrain_names = eigenstrainpv2
    block = 0
  [../]
  [./strain_phasegp_g1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasegp
    eigenstrain_names = eigenstrainpv2
    block = 1
  [../]
  [./strain_phasegpp1_g0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasegpp1
    eigenstrain_names = eigenstrainpv1
    block = 0
  [../]
  [./strain_phasegpp1_g1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasegpp1
    eigenstrain_names = eigenstrainpv1
    block = 1
  [../]
  
  [./strain_phasegpp2_g0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasegpp2
    eigenstrain_names = eigenstrainpv3
    block = 0
  [../]
  [./strain_phasegpp2_g1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasegpp2
    eigenstrain_names = eigenstrainpv3
    block = 1
  [../]
  [./strain_phased_g0]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phased
    eigenstrain_names = eigenstrainpv4
    block = 0
  [../]
  [./strain_phased_g1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phased
    eigenstrain_names = eigenstrainpv4
    block = 1
  [../]


[./eigen_straingp_g0]
    type = ComputeRotatedEigenstrain
    base_name = phasegp
    eigen_base = '-0.003 -0.003 0 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv2
    block = 0
  [../]
   [./eigen_straingp_g1]
    type = ComputeRotatedEigenstrain
    base_name = phasegp
    eigen_base = '-0.003 -0.003 0 0 0 0'
    Euler_angles = '45 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv2
    block = 1
  [../]
  [./eigen_strainpv1_g0]
    type = ComputeRotatedEigenstrain
    base_name = phasegpp1
    eigen_base = '0.028 0.0067 0 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv1
    block = 0
  [../]
  [./eigen_strainpv1_g1]
    type = ComputeRotatedEigenstrain
    base_name = phasegpp1
    eigen_base = '0.028 0.0067 0 0 0 0'
    Euler_angles = '45 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv1
    block = 1
  [../]
  [./eigen_strainpv3_g0]
    type = ComputeRotatedEigenstrain
    base_name = phasegpp2
    eigen_base = '0.0067 0.028 0 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv3
    block = 0
  [../]
  [./eigen_strainpv3_g1]
    type = ComputeRotatedEigenstrain
    base_name = phasegpp2
    eigen_base = '0.0067 0.028 0 0 0 0'
    Euler_angles = '45 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv3
    block = 1
  [../]
  # per-block eigenstrain (mirrors stiffness blocks)
[./eigen_strainpv4_g0]
  type = ComputeRotatedEigenstrain
  base_name = phased
  eigen_base = '0.005320 0.01332 0.02378 0 0 0'
  Euler_angles = '5.176036589 1.150261992 2.456873451'
  prefactor = 0
  eigenstrain_name = eigenstrainpv4
  block = 0
[../]

[./eigen_strainpv4_g1]
  type = ComputeRotatedEigenstrain
  base_name = phased
  eigen_base = '0.005320 0.01332 0.02378 0 0 0'
  Euler_angles = '5.176036589 1.150261992 2.456873451'
  prefactor = 0
  eigenstrain_name = eigenstrainpv4
  block = 1
[../]
  


  # Generate the global stress from the phase stresses
  [./global_stress]
    type = MultiPhaseStressMaterial
    phase_base = 'phasegp phasegpp1 phasegpp2 phased phasem'
    h          = 'hgp     hgpp1   hgpp2   hd   hm'
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
  
  # Kernels for Allen-Cahn equation for eta_gp
  [./deta_gp_dt]
    type = TimeDerivative
    variable = eta_gp
  [../]
  [./ACBulkF_gp]
    type = KKSMultiACBulkF
    variable  = eta_gp
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hd hm'
    gi_name   = gpv1
    eta_i     = eta_gp
    wi        = 0.01
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1d c1m c2gp c2gpp1 c2gpp2 c2d c2m eta_gpp1 eta_gpp2 eta_d eta_m'
  [../]
  [./ACBulkCpv1_c1]
    type = KKSMultiACBulkC
    variable  = eta_gp
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hd hm'
    cj_names  = 'c1gp c1gpp1 c1gpp2 c1d c1m'
    eta_i     = eta_gp
    coupled_variables      = 'c2gp c2gpp1 c2gpp2 c2d c2m eta_gpp1 eta_gpp2 eta_d eta_m'
  [../]
  [./ACBulkCpv1_c2]
    type = KKSMultiACBulkC
    variable  = eta_gp
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hd hm'
    cj_names  = 'c2gp c2gpp1 c2gpp2 c2d c2m'
    eta_i     = eta_gp
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1d c1m  eta_gpp1 eta_gpp2 eta_d eta_m'
  [../]
  [./ACInterfacepv1]
    type = ACInterface
    variable = eta_gp
    kappa_name = kappa
  [../]

  # Kernels for Allen-Cahn equation for eta_gpp1
  [./deta_gpp1_dt]
    type = TimeDerivative
    variable = eta_gpp1
  [../]
  [./ACBulkF_gpp1]
    type = KKSMultiACBulkF
    variable  = eta_gpp1
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hd hm'
    gi_name   = gpv2
    eta_i     = eta_gpp1
    wi        = 0.01
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1d c1m c2gp c2gpp1 c2gpp2 c2d c2m eta_gp eta_gpp2 eta_d eta_m'
  [../]
  [./ACBulkCpv2_c1]
    type = KKSMultiACBulkC
    variable  = eta_gpp1
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hd hm'
    cj_names  = 'c1gp c1gpp1 c1gpp2 c1d c1m'
    eta_i     = eta_gpp1
    coupled_variables      = 'c2gp c2gpp1 c2gpp2 c2d c2m eta_gp eta_gpp2 eta_d eta_m'
  [../]
  [./ACBulkCpv2_c2]
    type = KKSMultiACBulkC
    variable  = eta_gpp1
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hd hm'
    cj_names  = 'c2gp c2gpp1 c2gpp2 c2d c2m'
    eta_i     = eta_gpp1
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1d c1m eta_gp eta_gpp2 eta_d eta_m'
  [../]
  [./ACInterfacepv2]
    type = ACInterface
    variable = eta_gpp1
    kappa_name = kappa
  [../]

  # Kernels for Allen-Cahn equation for eta_gpp2
  [./deta_gpp2_dt]
    type = TimeDerivative
    variable = eta_gpp2
  [../]
  [./ACBulkF_gpp2]
    type = KKSMultiACBulkF
    variable  = eta_gpp2
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hd hm'
    gi_name   = gpv3
    eta_i     = eta_gpp2
    wi        = 0.01
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1d c1m c2gp c2gpp1 c2gpp2 c2d c2m eta_gp eta_gpp1 eta_d eta_m'
  [../]
  [./ACBulkCpv3_c1]
    type = KKSMultiACBulkC
    variable  = eta_gpp2
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hd hm'
    cj_names  = 'c1gp c1gpp1 c1gpp2 c1d c1m'
    eta_i     = eta_gpp2
    coupled_variables      = 'c2gp c2gpp1 c2gpp2 c2d c2m eta_gp eta_gpp1 eta_d eta_m'
  [../]
  [./ACBulkCpv3_c2]
    type = KKSMultiACBulkC
    variable  = eta_gpp2
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hd hm'
    cj_names  = 'c2gp c2gpp1 c2gpp2 c2d c2m'
    eta_i     = eta_gpp2
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1d c1m eta_gp eta_gpp1 eta_d eta_m'
  [../]
  [./ACInterfacepv3]
    type = ACInterface
    variable = eta_gpp2
    kappa_name = kappa
  [../]
 
# Kernels for Allen-Cahn equation for eta_d
  [./deta_d_dt]
    type = TimeDerivative
    variable = eta_d
  [../]
  [./ACBulkF_d]
    type = KKSMultiACBulkF
    variable  = eta_d
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hd hm'
    gi_name   = gpv4
    eta_i     = eta_d
    wi        = 0.01
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1d c1m c2gp c2gpp1 c2gpp2 c2d c2m eta_gp eta_gpp1 eta_gpp2 eta_m'
  [../]
  [./ACBulkCpv4_c1]
    type = KKSMultiACBulkC
    variable  = eta_d
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hd hm'
    cj_names  = 'c1gp c1gpp1 c1gpp2 c1d c1m'
    eta_i     = eta_d
    coupled_variables      = 'c2gp c2gpp1 c2gpp2 c2d c2m eta_gp eta_gpp1 eta_gpp2 eta_m'
  [../]
  [./ACBulkCpv4_c2]
    type = KKSMultiACBulkC
    variable  = eta_d
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hd hm'
    cj_names  = 'c2gp c2gpp1 c2gpp2 c2d c2m'
    eta_i     = eta_d
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1d c1m eta_gp eta_gpp1 eta_gpp2 eta_m'
  [../]
  [./ACInterfacepv4]
    type = ACInterface
    variable = eta_d
    kappa_name = kappa
  [../]

# Kernels for constraint equation |eta_gp| + |eta_gpp1| + eta_m = 1
  # eta3 is the nonlinear variable for the constraint equation
  [./eta_mreaction]
    type = MatReaction
    variable = eta_m
    reaction_rate = 1
 [../]
  [./eta_gpreaction]
    type = MatReaction_abscouple
    variable = eta_m
    v = eta_gp
    reaction_rate = 1
  [../]
  [./eta_gpp1reaction]
    type = MatReaction_abscouple
    variable = eta_m
    v = eta_gpp1
    reaction_rate = 1
  [../]
  [./eta_gpp2reaction]
    type = MatReaction_abscouple
    variable = eta_m
    v = eta_gpp2
    reaction_rate = 1
  [../]
  [./eta_dreaction]
    type = MatReaction_abscouple
    variable = eta_m
    v = eta_d
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
  [./diff_c1gp]
    type = MatDiffusion
    variable = c1
    diffusivity = Dhgp
    v = c1gp
  [../]
  [./diff_c1gpp1]
    type = MatDiffusion
    variable = c1
    diffusivity = Dhgpp1
    v = c1gpp1
  [../]
  [./diff_c1gpp2]
    type = MatDiffusion
    variable = c1
    diffusivity = Dhgpp2
    v = c1gpp2
  [../]
  [./diff_c1d]
    type = MatDiffusion
    variable = c1
    diffusivity = Dhd
    v = c1d
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
  [./diff_c2gp]
    type = MatDiffusion
    variable = c2
    diffusivity = Dhgp
    v = c2gp
  [../]
  [./diff_c2gpp1]
    type = MatDiffusion
    variable = c2
    diffusivity = Dhgpp1
    v = c2gpp1
  [../]
  [./diff_c2gpp2]
    type = MatDiffusion
    variable = c2
    diffusivity = Dhgpp2
    v = c2gpp2
  [../]   
  [./diff_c2d]
    type = MatDiffusion
    variable = c2
    diffusivity = Dhd
    v = c2d
  [../]

  # Phase concentration constraints
   [./chempot1m_pv1]
    type = KKSPhaseChemicalPotential
    variable = c1m
    cb       = c1gp
    fa_name  = F_gamma
    fb_name  = F_gp
    args_a   = c2m
    args_b   = c2gp
  [../]
 [./chempot1m_pv2]
    type = KKSPhaseChemicalPotential
    variable = c1gp
    cb       = c1gpp1
    fa_name  = F_gp
    fb_name  = F_gpp1
    args_a   = c2gp
    args_b   = c2gpp1
  [../]
  [./chempot1m_pv3]
    type = KKSPhaseChemicalPotential
    variable = c1gpp1
    cb       = c1gpp2
    fa_name  = F_gpp1
    fb_name  = F_gpp2
    args_a   = c2gpp1
    args_b   = c2gpp2
  [../]
  [./chempot1m_pv4]
    type = KKSPhaseChemicalPotential
    variable = c1gpp2
    cb       = c1d
    fa_name  = F_gpp2
    fb_name  = F_d
    args_a   = c2gpp2
    args_b   = c2d
  [../]
  [./chempot2m_pv1]
    type = KKSPhaseChemicalPotential
    variable = c2m
    cb       = c2gp
    fa_name  = F_gamma
    fb_name  = F_gp
    args_a   = c1m
    args_b   = c1gp
  [../]
  [./chempot2m_pv2]
    type = KKSPhaseChemicalPotential
    variable = c2gp
    cb       = c2gpp1
    fa_name  = F_gp
    fb_name  = F_gpp1
    args_a   = c1gp
    args_b   = c1gpp1
  [../]
  [./chempot2m_pv3]
    type = KKSPhaseChemicalPotential
    variable = c2gpp1
    cb       = c2gpp2
    fa_name  = F_gpp1
    fb_name  = F_gpp2
    args_a   = c1gpp1
    args_b   = c1gpp2
  [../]
  [./chempot2m_pv4]
    type = KKSPhaseChemicalPotential
    variable = c2gpp2
    cb       = c2d
    fa_name  = F_gpp2
    fb_name  = F_d
    args_a   = c1gpp2
    args_b   = c1d
  [../]
    
  [./phaseconcentration_c1gpp2]
    type = KKSMultiPhaseConcentration
    variable = c1gpp2
    cj = 'c1m c1gp c1gpp1 c1gpp2 c1d'
    hj_names = 'hm hgp hgpp1 hgpp2 hd'
    etas = 'eta_m eta_gp eta_gpp1 eta_gpp2 eta_d'
    c = c1
  [../]

  [./phaseconcentration_c1d]
    type = KKSMultiPhaseConcentration
    variable = c1d
    cj = 'c1m c1gp c1gpp1 c1gpp2 c1d'
    hj_names = 'hm hgp hgpp1 hgpp2 hd'
    etas = 'eta_m eta_gp eta_gpp1 eta_gpp2 eta_d'
    c = c1
  [../]

  [./phaseconcentration_c2gpp2]
    type = KKSMultiPhaseConcentration
    variable = c2gpp2
    cj = 'c2m c2gp c2gpp1 c2gpp2 c2d'
    hj_names = 'hm hgp hgpp1 hgpp2 hd'
    etas = 'eta_m eta_gp eta_gpp1 eta_gpp2 eta_d'
    c = c2
  [../]
    
  [./phaseconcentration_c2d]
    type = KKSMultiPhaseConcentration
    variable = c2d
    cj = 'c2m c2gp c2gpp1 c2gpp2 c2d'
    hj_names = 'hm hgp hgpp1 hgpp2 hd'
    etas = 'eta_m eta_gp eta_gpp1 eta_gpp2 eta_d'
    c = c2
  [../]
  

[]

[AuxKernels]


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
    Fj_names = 'F_gp F_gpp1 F_gpp2 F_d F_gamma'
    hj_names = 'hgp hgpp1 hgpp2 hd hm'
    gj_names = 'gpv1 gpv2 gpv3 gpv4 gm'
    variable = Energy
    w = 1
    interfacial_vars =  'eta_gp  eta_gpp1  eta_gpp2  eta_d  eta_m'
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
  property = hgp 
  execute_on = timestep_end 
  [../]
  [./copy_h_pv2] 
  type = MaterialRealAux 
  variable = h_pv2_aux 
  property = hgpp1 
  execute_on = timestep_end 
  [../]
  [./copy_h_pv3] 
  type = MaterialRealAux 
  variable = h_pv3_aux 
  property = hgpp2 
  execute_on = timestep_end 
  [../]
  [./copy_h_pv4]
  type = MaterialRealAux
  variable = h_pv4_aux
  property = hd
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
     variable = eta_gp
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
  
 # Average hpv in whole domain
    [./af_pv1] 
      type = ElementAverageMaterialProperty 
      mat_prop = hgp 
    [../]
    [./af_pv2] 
      type = ElementAverageMaterialProperty 
      mat_prop = hgpp1 
    [../]
    [./af_pv3] 
      type = ElementAverageMaterialProperty 
      mat_prop = hgpp2 
    [../]
  [./af_pv4]
    type = ElementAverageMaterialProperty
    mat_prop = hd
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
  [./eta_gp]
    type = ElementIntegralVariablePostprocessor
    variable = eta_gp
    use_absolute_value = true
  [../]
   [./eta_gpp1]
    type = ElementIntegralVariablePostprocessor
    variable = eta_gpp1
    use_absolute_value = true
  [../]
  [./eta_gpp2]
    type = ElementIntegralVariablePostprocessor
    variable = eta_gpp2
    use_absolute_value = true
  [../]
  [./eta_d]
    type = ElementIntegralVariablePostprocessor
    variable = eta_d
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
