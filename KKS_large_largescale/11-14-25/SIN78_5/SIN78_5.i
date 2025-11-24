# This test is for the multicomponent In718 alloy

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 175
  ny = 175
#  nz = 2
  xmin = 0
  xmax = 350
  ymin = 0
  ymax = 350
  zmin = 0
  zmax = 0
#  elem_type = QUAD4
[]
[BCs]
  [./u_right_pull]
    type = Pressure
    displacements = 'disp_x disp_y'
    variable = disp_x
    boundary = right
    factor = -0.2
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
  [./hgp_aux]
     family = MONOMIAL 
     order = CONSTANT 
  [../]
  [./hgpp1_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./hgpp2_aux]
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
  [./vonmises_hgp]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_hgpp1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_hgpp2]
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
    # New clamped auxiliary variables
  [./eta_gp_clamped]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./eta_gpp1_clamped]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./eta_gpp2_clamped]
    family = MONOMIAL
    order = CONSTANT
  [../]
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
# Normalized material constants (paper normalization applied)
# ------------------------------------------------------------------
# Normalization factor used in paper: 1.49e8 J/mol
# We set Vm_norm = 1.0 so expressions are already in paper-normalized units.


# ------------------------------------------------------------------
# 1. Chemical free energy of matrix (γ) — normalized
# ------------------------------------------------------------------
[./fc_gamma]
  type = DerivativeParsedMaterial
  property_name = fc_gamma
  coupled_variables = 'c1m c2m'               # c1m = Al, c2m = Nb (γ-phase)
  constant_names = 'Vm_norm C_Al_norm C_Nb_norm c1_eq_gamma c2_eq_gamma'
  constant_expressions = '1.0 50 50 0.0161 0.00723'
  expression = '( C_Al_norm*(c1m - c1_eq_gamma)^2 + C_Nb_norm*(c2m - c2_eq_gamma)^2 ) / Vm_norm'
[../]

[./elastic_free_energy_m]
  type = ElasticEnergyMaterial
  base_name = phasem
  property_name = fe_gamma
  coupled_variables = ' '
[../]

[./Total_energy_m]
  type = DerivativeSumMaterial
  property_name = F_gamma
  sum_materials = 'fc_gamma fe_gamma'
  coupled_variables = 'c1m c2m'
[../]

# ------------------------------------------------------------------
# 2. Chemical free energy of γ' (L12) — normalized
# ------------------------------------------------------------------
[./fc_gp]
  type = DerivativeParsedMaterial
  property_name = fc_gp
  coupled_variables = 'c1gp c2gp'
  constant_names = 'Vm_norm C_Al_norm C_Nb_norm c1_eq_gp c2_eq_gp'
  constant_expressions = '1.0 50 50 0.187 0.0157'
  expression = '( C_Al_norm*(c1gp - c1_eq_gp)^2 + C_Nb_norm*(c2gp - c2_eq_gp)^2 ) / Vm_norm'
[../]

[./elastic_free_energy_gp]
  type = ElasticEnergyMaterial
  base_name = phasegp
  property_name = fe_gp
  coupled_variables = ' '
[../]

[./Total_energy_gp]
  type = DerivativeSumMaterial
  property_name = F_gp
  sum_materials = 'fc_gp fe_gp'
  coupled_variables = 'c1gp c2gp'
[../]

# ------------------------------------------------------------------
# 3. γ'' variant 1 (metastable, "delta") — normalized (Δf included)
# ------------------------------------------------------------------
[./fc_gpp1]
  type = DerivativeParsedMaterial
  property_name = fc_gpp1
  coupled_variables = 'c1gpp1 c2gpp1'
  # Use same normalized C constants and the normalized Delta_f
  constant_names = 'Vm_norm C_Al_norm C_Nb_norm c1_eq_d c2_eq_d Delta_f_norm'
  constant_expressions = '1.0 50 50 7.27e-4 0.196 1.0392617e-5'
  expression = '( C_Al_norm*(c1gpp1 - c1_eq_d)^2 + C_Nb_norm*(c2gpp1 - c2_eq_d)^2 + Delta_f_norm ) / Vm_norm'
[../]

[./elastic_free_energy_gpp1]
  type = ElasticEnergyMaterial
  base_name = phasegpp1
  property_name = fe_gpp1
  coupled_variables = ' '
[../]

[./Total_energy_gpp1]
  type = DerivativeSumMaterial
  property_name = F_gpp1
  sum_materials = 'fc_gpp1 fe_gpp1'
  coupled_variables = 'c1gpp1 c2gpp1'
[../]

# ------------------------------------------------------------------
# 4. γ'' variant 2 (metastable, second delta variant) — normalized
# ------------------------------------------------------------------
[./fc_gpp2]
  type = DerivativeParsedMaterial
  property_name = fc_gpp2
  coupled_variables = 'c1gpp2 c2gpp2'
  constant_names = 'Vm_norm C_Al_norm C_Nb_norm c1_eq_d c2_eq_d Delta_f_norm'
  constant_expressions = '1.0 50 50 7.27e-4 0.196 1.0392617e-5'
  expression = '( C_Al_norm*(c1gpp2 - c1_eq_d)^2 + C_Nb_norm*(c2gpp2 - c2_eq_d)^2 + Delta_f_norm ) / Vm_norm'
[../]

[./elastic_free_energy_gpp2]
  type = ElasticEnergyMaterial
  base_name = phasegpp2
  property_name = fe_gpp2
  coupled_variables = ' '
[../]

[./Total_energy_gpp2]
  type = DerivativeSumMaterial
  property_name = F_gpp2
  sum_materials = 'fc_gpp2 fe_gpp2'
  coupled_variables = 'c1gpp2 c2gpp2'
[../]


    # Switching functions for each phase
  # hm(eta_gp, eta_gpp1, eta_m)
  [./hm]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_m
    all_etas = 'eta_gp eta_gpp1 eta_gpp2 eta_m'
    h_name = hm
  [../]
  # hgp(eta_gp, eta_gpp1, eta_m)
  [./hgp]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_gp
    all_etas = 'eta_gp eta_gpp1 eta_gpp2 eta_m'
    h_name = hgp
  [../]
  # hgpp1(eta_gp, eta_gpp1, eta_m)
  [./hgpp1]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_gpp1
    all_etas = 'eta_gp eta_gpp1 eta_gpp2 eta_m'
    h_name = hgpp1
  [../]
[./hgpp2]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_gpp2
    all_etas = 'eta_gp eta_gpp1 eta_gpp2 eta_m'
    h_name = hgpp2
[../]
  # Coefficients for diffusion equation
  [./Dhm]
    type = DerivativeParsedMaterial
    material_property_names = 'D hm'
    expression = D*hm
    property_name = Dhm
  [../]
  [./Dhgp]
    type = DerivativeParsedMaterial
    material_property_names = 'D hgp'
    expression = D*hgp
    property_name = Dhgp
  [../]
  [./Dhgpp1]
    type = DerivativeParsedMaterial
    material_property_names = 'D hgpp1'
    expression = D*hgpp1
    property_name = Dhgpp1
  [../]
  [./Dhgpp2]
    type = DerivativeParsedMaterial
    material_property_names = 'D hgpp2'
    expression = D*hgpp2
    property_name = Dhgpp2
  [../]

# Barrier functions for each phase
  [./gm]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta_m
    function_name = gm
  [../]
  [./gp]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_gp
    function_name = gp
  [../]
  [./gpp1]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_gpp1
    function_name = gpp1
  [../]
  [./gpp2]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta_gpp2
    function_name = gpp2
  [../]

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'L    kappa  D   misfit     W'
    prop_values = '0.3  0.01   1     1      0.01'
  [../]
 
#Mechanical properties
  [./Stiffness_phasem]
    type = ComputeElasticityTensor
    C_ijkl = '272.1 169 169 272.1 169 272.1 131 131 131' #Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of     
    base_name = phasem
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    fill_method = symmetric9
  [../]
   [./Stiffness_phasegp]
    type = ComputeElasticityTensor
    C_ijkl = '243 154.8 154.8 243 154.8 243 132.3 132.3 132.3'
    base_name = phasegp
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    fill_method = symmetric9
  [../]
  [./Stiffness_phasegpp1]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'#Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of    
    base_name = phasegpp1
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    fill_method = symmetric9
  [../]
  [./Stiffness_phasegpp2]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
    base_name = phasegpp2
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    fill_method = symmetric9
  [../]

  [./stress_phasegp]
    type = ComputeLinearElasticStress
    base_name = phasegp
  [../]
  [./stress_phasegpp1]
    type = ComputeLinearElasticStress
    base_name = phasegpp1
  [../]
  [./stress_phasegpp2]
    type = ComputeLinearElasticStress
    base_name = phasegpp2
  [../]
  [./stress_phasem]
    type = ComputeLinearElasticStress
    base_name = phasem
  [../]

  [./strain_phasem]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasem
  [../]
  [./strain_phasegp]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasegp
    eigenstrain_names = eigenstrainpv1
  [../]
  [./strain_phasegpp1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasegpp1
    eigenstrain_names = 'eigenstrainpv2'
  [../]
  [./strain_phasegpp2]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasegpp2
    eigenstrain_names = 'eigenstrainpv3'
  [../]
  [./eigen_strainpv1]
    type = ComputeRotatedEigenstrain
    base_name = phasegp
    eigen_base = '-0.003 -0.003 -0.003 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv1
  [../]
 [./eigen_strainpv2]
    type = ComputeRotatedEigenstrain
    base_name = phasegpp1
    eigen_base = '0.0286 0.0067 0.0067 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv2
  [../]
 [./eigen_strainpv3]
    type = ComputeRotatedEigenstrain
    base_name = phasegpp2
    eigen_base = '0.0067 0.0286 0.0067 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv3
  [../]


  # Generate the global stress from the phase stresses
  [./global_stress]
    type = MultiPhaseStressMaterial
    phase_base = 'phasegp phasegpp1 phasegpp2 phasem'
    h          = 'hgp     hgpp1   hgpp2   hm'
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
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hm'
    gi_name   = gp
    eta_i     = eta_gp
    wi        = 0.01
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1m c2gp c2gpp1 c2gpp2 c2m eta_gpp1 eta_gpp2 eta_m'
  [../]
  [./ACBulkCpv1_c1]
    type = KKSMultiACBulkC
    variable  = eta_gp
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hm'
    cj_names  = 'c1gp c1gpp1 c1gpp2 c1m'
    eta_i     = eta_gp
    coupled_variables      = 'c2gp c2gpp1 c2gpp2 c2m eta_gpp1 eta_gpp2 eta_m'
  [../]
  [./ACBulkCpv1_c2]
    type = KKSMultiACBulkC
    variable  = eta_gp
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hm'
    cj_names  = 'c2gp c2gpp1 c2gpp2 c2m'
    eta_i     = eta_gp
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1m eta_gpp1 eta_gpp2 eta_m'
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
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hm'
    gi_name   = gpp1
    eta_i     = eta_gpp1
    wi        = 0.01
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1m c2gp c2gpp1 c2gpp2 c2m eta_gp eta_gpp2 eta_m'
  [../]
  [./ACBulkCpv2_c1]
    type = KKSMultiACBulkC
    variable  = eta_gpp1
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hm'
    cj_names  = 'c1gp c1gpp1 c1gpp2 c1m'
    eta_i     = eta_gpp1
    coupled_variables      = 'c2gp c2gpp1 c2gpp2 c2m eta_gp eta_gpp2 eta_m'
  [../]
  [./ACBulkCpv2_c2]
    type = KKSMultiACBulkC
    variable  = eta_gpp1
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hm'
    cj_names  = 'c2gp c2gpp1 c2gpp2 c2m'
    eta_i     = eta_gpp1
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1m eta_gp eta_gpp2 eta_m'
  [../]
  [./ACInterfacepv2]
    type = ACInterface
    variable = eta_gpp1
    kappa_name = kappa
  [../]

  # Kernels for Allen-Cahn equation for eta_gpp1
  [./deta_gpp2_dt]
    type = TimeDerivative
    variable = eta_gpp2
  [../]
  [./ACBulkF_gpp2]
    type = KKSMultiACBulkF
    variable  = eta_gpp2
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hm'
    gi_name   = gpp2
    eta_i     = eta_gpp2
    wi        = 0.01
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1m c2gp c2gpp1 c2gpp2 c2m eta_gp eta_gpp1 eta_m'
  [../]
  [./ACBulkCpv3_c1]
    type = KKSMultiACBulkC
    variable  = eta_gpp2
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hm'
    cj_names  = 'c1gp c1gpp1 c1gpp2 c1m'
    eta_i     = eta_gpp2
    coupled_variables      = 'c2gp c2gpp1 c2gpp2 c2m eta_gp eta_gpp1 eta_m'
  [../]
  [./ACBulkCpv3_c2]
    type = KKSMultiACBulkC
    variable  = eta_gpp2
    Fj_names  = 'F_gp F_gpp1 F_gpp2 F_gamma'
    hj_names  = 'hgp hgpp1 hgpp2 hm'
    cj_names  = 'c2gp c2gpp1 c2gpp2 c2m'
    eta_i     = eta_gpp2
    coupled_variables      = 'c1gp c1gpp1 c1gpp2 c1m eta_gp eta_gpp1 eta_m'
  [../]
  [./ACInterfacepv3]
    type = ACInterface
    variable = eta_gpp2
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
  [./one]
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
   
  [./phaseconcentration_c1m]
    type = KKSMultiPhaseConcentration
    variable = c1m
    cj = 'c1m c1gp c1gpp1 c1gpp2'
    hj_names = 'hm hgp hgpp1 hgpp2'
    etas = 'eta_m eta_gp eta_gpp1 eta_gpp2'
    c = c1
  [../]
  [./phaseconcentration_c2m]
    type = KKSMultiPhaseConcentration
    variable = c2m
    cj = 'c2m c2gp c2gpp1 c2gpp2'
    hj_names = 'hm hgp hgpp1 hgpp2'
    etas = 'eta_m eta_gp eta_gpp1 eta_gpp2'
    c = c2
  [../]
  [./phaseconcentration_c1gp]
    type = KKSMultiPhaseConcentration
    variable = c1gp
    cj = 'c1m c1gp c1gpp1 c1gpp2'
    hj_names = 'hm hgp hgpp1 hgpp2'
    etas = 'eta_m eta_gp eta_gpp1 eta_gpp2'
    c = c1
  [../]
  [./phaseconcentration_c2gp]
    type = KKSMultiPhaseConcentration
    variable = c2gp
    cj = 'c2m c2gp c2gpp1 c2gpp2'
    hj_names = 'hm hgp hgpp1 hgpp2'
    etas = 'eta_m eta_gp eta_gpp1 eta_gpp2'
    c = c2
  [../]
  [./phaseconcentration_c1gpp1]
    type = KKSMultiPhaseConcentration
    variable = c1gpp1
    cj = 'c1m c1gp c1gpp1 c1gpp2'
    hj_names = 'hm hgp hgpp1 hgpp2'
    etas = 'eta_m eta_gp eta_gpp1 eta_gpp2'
    c = c1
  [../]
  [./phaseconcentration_c2gpp1]
    type = KKSMultiPhaseConcentration
    variable = c2gpp1
    cj = 'c2m c2gp c2gpp1 c2gpp2'
    hj_names = 'hm hgp hgpp1 hgpp2'
    etas = 'eta_m eta_gp eta_gpp1 eta_gpp2'
    c = c2
  [../]
  [./phaseconcentration_c1gpp2]
    type = KKSMultiPhaseConcentration
    variable = c1gpp2
    cj = 'c1m c1gp c1gpp1 c1gpp2'
    hj_names = 'hm hgp hgpp1 hgpp2'
    etas = 'eta_m eta_gp eta_gpp1 eta_gpp2'
    c = c1
  [../]
  [./phaseconcentration_c2gpp2]
    type = KKSMultiPhaseConcentration
    variable = c2gpp2
    cj = 'c2m c2gp c2gpp1 c2gpp2'
    hj_names = 'hm hgp hgpp1 hgpp2'
    etas = 'eta_m eta_gp eta_gpp1 eta_gpp2'
    c = c2
  [../]
  

[]

[AuxKernels]
  [temperature]
    type = FunctionAux
    variable = temperature
    function = '1073'
    execute_on = timestep_begin
  []
  [./Energy_total]
    type = KKSMultiFreeEnergy
    Fj_names = 'F_gp F_gpp1 F_gpp2 F_gamma'
    hj_names = 'hgp hgpp1 hgpp2 hm'
    gj_names = 'gp gpp1 gpp2 gm'
    variable = Energy
    w = 1
    interfacial_vars =  'eta_gp  eta_gpp1  eta_gpp2  eta_m'
    kappa_names =       'kappa kappa kappa kappa'
  [../]
  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
    execute_on = timestep_end
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
  variable = hgp_aux 
  property = hgp 
  execute_on = timestep_end 
  [../]
  [./copy_h_pv2] 
  type = MaterialRealAux 
  variable = hgpp1_aux 
  property = hgpp1 
  execute_on = timestep_end 
  [../]
  [./copy_h_pv3] 
  type = MaterialRealAux 
  variable = hgpp2_aux 
  property = hgpp2 
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
     variable = vonmises_hgp
     coupled_variables = 'vonmises hgp_aux'
     expression = 'vonmises * hgp_aux'
     execute_on = timestep_end
   [../]
  [./vonmises_times_h_pv2]
      type = ParsedAux
      variable = vonmises_hgpp1
      coupled_variables = 'vonmises hgpp1_aux'
      expression = 'vonmises * hgpp1_aux'
      execute_on = timestep_end
    [../]
    [./vonmises_times_h_pv3]
      type = ParsedAux
      variable = vonmises_hgpp2
      coupled_variables = 'vonmises hgpp2_aux'
      expression = 'vonmises * hgpp2_aux'
      execute_on = timestep_end
    [../]
    [./vonmises_times_h_m]
      type = ParsedAux
      variable = vonmises_h_m
      coupled_variables = 'vonmises h_m_aux'
      expression = 'vonmises * h_m_aux'
      execute_on = timestep_end
    [../]
  
  # New clamped auxiliary kernels
  [./clamp_eta_gp]
    type = ParsedAux
    variable = eta_gp_clamped
    coupled_variables = 'eta_gp'
    expression = 'max(-1, min(1, eta_gp))'
    execute_on = timestep_end
  [../]
  [./clamp_eta_gpp1]
    type = ParsedAux
    variable = eta_gpp1_clamped
    coupled_variables = 'eta_gpp1'
    expression = 'max(-1, min(1, eta_gpp1))'
    execute_on = timestep_end
  [../]
  [./clamp_eta_gpp2]
    type = ParsedAux
    variable = eta_gpp2_clamped
    coupled_variables = 'eta_gpp2'
    expression = 'max(-1, min(1, eta_gpp2))'
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

  end_time = 100000

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
   [temperature]
    type = ElementAverageValue
    variable = temperature
  []
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
  variable = vonmises_hgp
  [../]
  [./num_vonmises_pv2]
  type = ElementIntegralVariablePostprocessor
  variable = vonmises_hgpp1
  [../]
  [./num_vonmises_pv3]
  type = ElementIntegralVariablePostprocessor
  variable = vonmises_hgpp2
  [../]
  [./num_vonmises_m]  
   type = ElementIntegralVariablePostprocessor
   variable = vonmises_h_m   
   [../]

  # Denominators: ∫ h dV  (phase volumes)
  [./den_pv1] 
  type = ElementIntegralVariablePostprocessor
  variable = hgp_aux
  use_absolute_value = true
  [../]
  [./den_pv2]
  type = ElementIntegralVariablePostprocessor
  variable = hgpp1_aux
  use_absolute_value = true
  [../]
  [./den_pv3] 
  type = ElementIntegralVariablePostprocessor 
  variable = hgpp2_aux 
  use_absolute_value = true
  [../]
  [./den_m]
   type = ElementIntegralVariablePostprocessor
   variable = h_m_aux
   use_absolute_value = true
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
    # Modified: Integrate clamped eta_pvX (instead of raw eta_pvX with absolute value)
  [./eta_gp_clamped]
    type = ElementIntegralVariablePostprocessor
    variable = eta_gp_clamped
    use_absolute_value = true
  [../]
  [./eta_gpp1_clamped]
    type = ElementIntegralVariablePostprocessor
    variable = eta_gpp1_clamped
    use_absolute_value = true 
  [../]
  [./eta_gpp2_clamped]
    type = ElementIntegralVariablePostprocessor
    variable = eta_gpp2_clamped
    use_absolute_value = true
  [../]
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
   

  # Time for plotting against the ramp
  [./time_pp]
    type = TimePostprocessor
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

