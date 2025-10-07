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
  [./vonmises_h_m]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx_h_pv1] 
    family = MONOMIAL 
    order = CONSTANT 
  [../]
  [./stress_xx_h_pv2]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./stress_xx_h_pv3]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./stress_xx_h_m]
    family = MONOMIAL
    order = CONSTANT
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
  [./eta_pv1]
    variable = eta_pv1
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 324
  [../]
  [./eta_pv2]
    variable = eta_pv2
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 230	
  [../]
  [./eta_pv3]
    variable = eta_pv3
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
  # Coefficients for diffusion equation
  [./Dhm]
    type = DerivativeParsedMaterial
    material_property_names = 'D hm'
    expression = D*hm
    property_name = Dhm
  [../]
  [./Dhpv1]
    type = DerivativeParsedMaterial
    material_property_names = 'D hpv1'
    expression = D*hpv1
    property_name = Dhpv1
  [../]
  [./Dhpv2]
    type = DerivativeParsedMaterial
    material_property_names = 'D hpv2'
    expression = D*hpv2
    property_name = Dhpv2
  [../]
  [./Dhpv3]
    type = DerivativeParsedMaterial
    material_property_names = 'D hpv3'
    expression = D*hpv3
    property_name = Dhpv3
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

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'L    kappa  D  misfit     W'
    prop_values = '0.3  0.01   1    1        0.01'
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
  [./Stiffness_phasepv1]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'#Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of    
    base_name = phasepv1
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    fill_method = symmetric9
  [../]
  [./Stiffness_phasepv2]
    type = ComputeElasticityTensor
    C_ijkl = '243 154.8 154.8 243 154.8 243 132.3 132.3 132.3'
    base_name = phasepv2
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    fill_method = symmetric9
  [../]
  [./Stiffness_phasepv3]
    type = ComputeElasticityTensor
    C_ijkl = '290.6 187 160.7 290.6 187 309.6 114.2 114.2 119.2'
    base_name = phasepv3
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    fill_method = symmetric9
  [../]

  [./stress_phasepv1]
    type = ComputeLinearElasticStress
    base_name = phasepv1
  [../]
  [./stress_phasepv2]
    type = ComputeLinearElasticStress
    base_name = phasepv2
  [../]
  [./stress_phasepv3]
    type = ComputeLinearElasticStress
    base_name = phasepv3
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
  [./strain_phasepv1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv1
    eigenstrain_names = eigenstrainpv1
  [../]
  [./strain_phasepv2]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv2
    eigenstrain_names = eigenstrainpv2
  [../]
  [./strain_phasepv3]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv3
    eigenstrain_names = eigenstrainpv3
  [../]

  [./eigen_strainpv1]
    type = ComputeRotatedEigenstrain
    base_name = phasepv1
    eigen_base = '0.01 0.001 0 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv1
  [../]

  [./eigen_strainpv2]
    type = ComputeRotatedEigenstrain
    base_name = phasepv2
    eigen_base = '-0.001 -0.001 0 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv2
  [../]
  [./eigen_strainpv3]
    type = ComputeRotatedEigenstrain
    base_name = phasepv3
    eigen_base = '0.001 0.01 0 0 0 0'
    Euler_angles = '0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv3
  [../]


  # Generate the global stress from the phase stresses
  [./global_stress]
    type = MultiPhaseStressMaterial
    phase_base = 'phasepv1 phasepv2 phasepv3 phasem'
    h          = 'hpv1     hpv2   hpv3   hm'
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
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hm'
    gi_name   = gpv1
    eta_i     = eta_pv1
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m c2pv1 c2pv2 c2pv3 c2m eta_pv2 eta_pv3 eta_m'
  [../]
  [./ACBulkCpv1_c1]
    type = KKSMultiACBulkC
    variable  = eta_pv1
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hm'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1m'
    eta_i     = eta_pv1
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2m eta_pv2 eta_pv3 eta_m'
  [../]
  [./ACBulkCpv1_c2]
    type = KKSMultiACBulkC
    variable  = eta_pv1
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hm'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2m'
    eta_i     = eta_pv1
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m eta_pv2 eta_pv3 eta_m'
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
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hm'
    gi_name   = gpv2
    eta_i     = eta_pv2
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m c2pv1 c2pv2 c2pv3 c2m eta_pv1 eta_pv3 eta_m'
  [../]
  [./ACBulkCpv2_c1]
    type = KKSMultiACBulkC
    variable  = eta_pv2
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hm'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1m'
    eta_i     = eta_pv2
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2m eta_pv1 eta_pv3 eta_m'
  [../]
  [./ACBulkCpv2_c2]
    type = KKSMultiACBulkC
    variable  = eta_pv2
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hm'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2m'
    eta_i     = eta_pv2
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m eta_pv1 eta_pv3 eta_m'
  [../]
  [./ACInterfacepv2]
    type = ACInterface
    variable = eta_pv2
    kappa_name = kappa
  [../]

  # Kernels for Allen-Cahn equation for eta_pv2
  [./deta_pv3_dt]
    type = TimeDerivative
    variable = eta_pv3
  [../]
  [./ACBulkFpv3]
    type = KKSMultiACBulkF
    variable  = eta_pv3
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hm'
    gi_name   = gpv3
    eta_i     = eta_pv3
    wi        = 0.01
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m c2pv1 c2pv2 c2pv3 c2m eta_pv1 eta_pv2 eta_m'
  [../]
  [./ACBulkCpv3_c1]
    type = KKSMultiACBulkC
    variable  = eta_pv3
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hm'
    cj_names  = 'c1pv1 c1pv2 c1pv3 c1m'
    eta_i     = eta_pv3
    coupled_variables      = 'c2pv1 c2pv2 c2pv3 c2m eta_pv1 eta_pv2 eta_m'
  [../]
  [./ACBulkCpv3_c2]
    type = KKSMultiACBulkC
    variable  = eta_pv3
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hm'
    cj_names  = 'c2pv1 c2pv2 c2pv3 c2m'
    eta_i     = eta_pv3
    coupled_variables      = 'c1pv1 c1pv2 c1pv3 c1m eta_pv1 eta_pv2 eta_m'
  [../]
  [./ACInterfacepv3]
    type = ACInterface
    variable = eta_pv3
    kappa_name = kappa
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
    
  [./phaseconcentration_c1pv1]
    type = KKSMultiPhaseConcentration
    variable = c1pv3
    cj = 'c1m c1pv1 c1pv2 c1pv3'
    hj_names = 'hm hpv1 hpv2 hpv3'
    etas = 'eta_m eta_pv1 eta_pv2 eta_pv3'
    c = c1
  [../]
  [./phaseconcentration_c2pv1]
    type = KKSMultiPhaseConcentration
    variable = c2pv3
    cj = 'c2m c2pv1 c2pv2 c2pv3'
    hj_names = 'hm hpv1 hpv2 hpv3'
    etas = 'eta_m eta_pv1 eta_pv2 eta_pv3'
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
    Fj_names = 'Fpv1 Fpv2 Fm'
    hj_names = 'hpv1 hpv2 hm'
    gj_names = 'gpv1 gpv2 gm'
    variable = Energy
    w = 1
    interfacial_vars =  'eta_pv1  eta_pv2  eta_m'
    kappa_names =       'kappa kappa kappa'
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
    [./vonmises_times_h_m]
      type = ParsedAux
      variable = vonmises_h_m
      coupled_variables = 'vonmises h_m_aux'
      expression = 'vonmises * h_m_aux'
      execute_on = timestep_end
    [../]
  [./sxx_times_h_pv1]
    type = ParsedAux
    variable = stress_xx_h_pv1
    coupled_variables = 'stress_xx h_pv1_aux'
    expression = 'stress_xx * h_pv1_aux'
    execute_on = timestep_end
  [../]
  [./sxx_times_h_pv2]
    type = ParsedAux
    variable = stress_xx_h_pv2
    coupled_variables = 'stress_xx h_pv2_aux'
    expression = 'stress_xx * h_pv2_aux'
    execute_on = timestep_end
  [../]
  [./sxx_times_h_pv3]
    type = ParsedAux
    variable = stress_xx_h_pv3
    coupled_variables = 'stress_xx h_pv3_aux'
    expression = 'stress_xx * h_pv3_aux'
    execute_on = timestep_end
  [../]
  [./sxx_times_h_m]
    type = ParsedAux
    variable = stress_xx_h_m
    coupled_variables = 'stress_xx h_m_aux'
    expression = 'stress_xx * h_m_aux'
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
  [./Adaptivity]
    initial_adaptivity = 0
    refine_fraction = 0.7
    coarsen_fraction = 0.1
    max_h_level = 1
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
  [./num_vonmises_m]  
   type = ElementIntegralVariablePostprocessor
   variable = vonmises_h_m   
   [../]
  [./num_pv1] 
  type = ElementIntegralVariablePostprocessor 
  variable = stress_xx_h_pv1 
  [../]
  [./num_pv2] 
  type = ElementIntegralVariablePostprocessor 
  variable = stress_xx_h_pv2 
  [../]
  [./num_pv3] 
  type = ElementIntegralVariablePostprocessor 
  variable = stress_xx_h_pv3 
  [../]
  [./num_m]  
   type = ElementIntegralVariablePostprocessor 
   variable = stress_xx_h_m   
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
  [./den_m]   
   type = ElementIntegralVariablePostprocessor 
   variable = h_m_aux   
   [../]
    
  [./eta_pv1]
    type = ElementIntegralVariablePostprocessor_new2
    variable = eta_pv1
  [../]
   [./eta_pv2]
    type = ElementIntegralVariablePostprocessor_new2
    variable = eta_pv2
  [../]
  [./eta_pv3]
    type = ElementIntegralVariablePostprocessor_new2
    variable = eta_pv3
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
