#
# KKS ternary (3 chemical component) system example in the split form
# We track c1 and c2 only, since c1 + c2 + c3 = 1
#

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 150
  ny = 150
#  nz = 2
  xmin = 0
  xmax = 300
  ymin = 0
  ymax = 300
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
  [./Energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./e_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Variables]
  # order parameter
  [./eta]
    order = FIRST
    family = LAGRANGE
  [../]
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



  # chemical potential solute 1
  [./w1]
    order = FIRST
    family = LAGRANGE
  [../]

  # chemical potential solute 2
  [./w2]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]


[]

[Functions]
  [./ic_func_eta]
    type = ParsedFunction
    expression = '0.5*(1.0-tanh((x)/sqrt(2.0)))'
  [../]
  [./ic_func_c1]
    type = ParsedFunction
    expression = '0.8*(0.5*(1.0-tanh(x/sqrt(2.0))))^3*(6*(0.5*(1.0-tanh(x/sqrt(2.0))))^2-15*(0.5*(1.0-tanh(x/sqrt(2.0))))+10)+0.1*(1-(0.5*(1.0-tanh(x/sqrt(2.0))))^3*(6*(0.5*(1.0-tanh(x/sqrt(2.0))))^2-15*(0.5*(1.0-tanh(x/sqrt(2.0))))+10))'
  [../]
  [./ic_func_c2]
    type = ParsedFunction
    expression = '0.1*(0.5*(1.0-tanh(x/sqrt(2.0))))^3*(6*(0.5*(1.0-tanh(x/sqrt(2.0))))^2-15*(0.5*(1.0-tanh(x/sqrt(2.0))))+10)+0.05*(1-(0.5*(1.0-tanh(x/sqrt(2.0))))^3*(6*(0.5*(1.0-tanh(x/sqrt(2.0))))^2-15*(0.5*(1.0-tanh(x/sqrt(2.0))))+10))'
  [../]

[]


[ICs]
  [./eta]
    variable = eta
    type = RandomIC
    min = -0.1625
    max = 0.1625
    seed = 192
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
  [./fm]
    type = DerivativeParsedMaterial
    f_name = fc_m
    args = 'c1m c2m'
    function = '50.0*((c1m-0.0161)^2+2*(c2m-0.00723)^2)'
  [../]
  # Elastic energy of the phase 0
  [./elastic_free_energy_m]
    type = ElasticEnergyMaterial
    base_name = phasem
    f_name = fe_m
    args = ' '
  [../]
  # Total free energy of the phase 0
  [./Total_energy_m]
    type = DerivativeSumMaterial
    f_name = Fm
    sum_materials = 'fc_m fe_m'
    args = 'c1m c2m'
  [../]

  [./fc_pv1]
    type = DerivativeParsedMaterial
    f_name = fc_pv1
    args = 'c1pv1 c2pv1'
    function = '50.0*((c1pv1-0.00727)^2+2*(c2pv1-0.196)^2)-2.0'
  [../]
  # Elastic energy of the phase 1
  [./elastic_free_energy_pv1]
    type = ElasticEnergyMaterial
    base_name = phasepv1
    f_name = fe_pv1
    args = ' '
  [../]
  # Total free energy of the phase 1
  [./Total_energy_pv1]
    type = DerivativeSumMaterial
    f_name = Fpv1
    sum_materials = 'fc_pv1 fe_pv1'
    args = 'c1pv1 c2pv1'
  [../]

  [./f2]
    type = DerivativeParsedMaterial
    f_name = fc_pv2
    args = 'c1pv2 c2pv2'
    function = '50.0*((c1pv2-0.187)^2+2*(c2pv2-0.000727)^2)'
  [../]
  # Elastic energy of the phase 2
  [./elastic_free_energy_pv2]
    type = ElasticEnergyMaterial
    base_name = phasepv2
    f_name = fe_pv2
    args = ' '
  [../]
  # Total free energy of the phase 2
  [./Total_energy_pv2]
    type = DerivativeSumMaterial
    f_name = Fpv2
    sum_materials = 'fc_pv2 fe_pv2'
    args = 'c1pv2 c2pv2'
  [../]

  # h(eta)
  [./h_eta]
    type = SwitchingFunctionMaterial
    h_order = HIGH
    eta = eta
  [../]

  # g(eta)
  [./g_eta]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta
  [../]

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'M   L   eps_sq kappa misfit'
    prop_values = '0.7 0.7 1.0    0.01 0000'
  [../]
  [./Stiffness_phasem]
    type = ComputeElasticityTensor
    C_ijkl = '2008.5 448.9 918.5 2008.5 918.5 1538.7 779.8 779.8 309.9'    
    base_name = phasem
    fill_method = symmetric9
  [../]
  [./Stiffness_phasepv1]
    type = ComputeElasticityTensor
    C_ijkl = '1982 534 496 1606 605 1788 985 1056 388'    
    base_name = phasepv1
    fill_method = symmetric9
  [../]
  [./Stiffness_phasepv2]
    type = ComputeElasticityTensor
    C_ijkl = '1606 534 605 1982 496 1788 1056 985 388'
    base_name = phasepv2
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
  
  [./eigen_strainpv1]
    type = ComputeEigenstrain
    base_name = phasepv1
    eigen_base = '1 0 0 0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv1
  [../]
  
  [./eigen_strainpv2]
    type = ComputeEigenstrain
    base_name = phasepv2
    eigen_base = '0 1 0 0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv2
  [../]
  
  
  # Generate the global stress from the phase stresses
  [./global_stress]
    type = MultiPhaseStressMaterial
    phase_base = 'phasepv1 phasepv2 phasem'
    h          = 'h     h     h'
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
  [./phaseconcentration_c1pv1]
    type = KKSMultiPhaseConcentration
    variable = c1pv1
    cj = 'c1m c1pv1 c1pv2'
    hj_names = 'h h h'
    etas = 'eta eta eta'
    c = c1
  [../]
  [./phaseconcentration_c2pv1]
    type = KKSMultiPhaseConcentration
    variable = c2pv1
    cj = 'c2m c2pv1 c2pv2'
    hj_names = 'h h h'
    etas = 'eta eta eta'
    c = c2
  [../]

  # Phase concentration constraints
  [./chempot1m_pv1]
    type = KKSPhaseChemicalPotential
    variable = c1m
    cb       = c1pv1
    fa_name  = Fm
    fb_name  = Fpv1
  [../]
  [./chempot1m_pv2]
    type = KKSPhaseChemicalPotential
    variable = c1pv2
    cb       = c1m
    fa_name  = Fpv2
    fb_name  = Fm
  [../]
  [./chempot2m_pv1]
    type = KKSPhaseChemicalPotential
    variable = c2m
    cb       = c2pv1
    fa_name  = Fm
    fb_name  = Fpv1
  [../]
  [./chempot2m_pv2]
    type = KKSPhaseChemicalPotential
    variable = c2pv2
    cb       = c2m
    fa_name  = Fpv2
    fb_name  = Fm
  [../]

 # Cahn-Hilliard Equations
  #
  [./CHBulk1_pv1]
    type = KKSSplitCHCRes
    variable = c1
    ca       = c1pv1
    fa_name  = fc_pv1
    w        = w1
    args_a   = 'c2pv1'
  [../]
  [./CHBulk2_pv1]
    type = KKSSplitCHCRes
    variable = c2
    ca       = c2pv1
    fa_name  = fc_pv1
    w        = w2
    args_a   = 'c1pv1'
  [../] 
  [./CHBulk1_pv2]
    type = KKSSplitCHCRes
    variable = c2
    ca       = c1pv2
    fa_name  = fc_pv2
    w        = w1
    args_a   = 'c2pv2'
  [../] 
  [./CHBulk2_pv2]
    type = KKSSplitCHCRes
    variable = c2
    ca       = c2pv2
    fa_name  = fc_pv2
    w        = w2
    args_a   = 'c1pv2'
  [../]  
  [./dc1dt]
    type = CoupledTimeDerivative
    variable = w1
    v = c1
  [../]
  [./dc2dt]
    type = CoupledTimeDerivative
    variable = w2
    v = c2
  [../]

  [./w1kernel]
    type = SplitCHWRes
    mob_name = M
    variable = w1
  [../]
  [./w2kernel]
    type = SplitCHWRes
    mob_name = M
    variable = w2
  [../]


  # Kernels for Allen-Cahn equation for eta_pv1
  [./deta_pv1_dt]
    type = TimeDerivative
    variable = eta
  [../]
  [./ACBulkFpv1]
    type = KKSMultiACBulkF
    variable  = eta
    Fj_names  = 'Fpv1 Fpv2 Fm'
    hj_names  = 'h h h'
    gi_name   = g
    eta_i     = eta
    wi        = 0.01
    args      = 'c1pv1 c1pv2 c1m c2pv1 c2pv2 c2m eta'
  [../]
  [./ACBulkCpv1_c1]
    type = KKSMultiACBulkC
    variable  = eta
    Fj_names  = 'Fpv1 Fpv2 Fm'
    hj_names  = 'h h h'
    cj_names  = 'c1pv1 c1pv2 c1m'
    eta_i     = eta
    args      = 'c2m c2pv1 c2pv2 eta'
  [../]
  [./ACBulkCpv1_c2]
    type = KKSMultiACBulkC
    variable  = eta
    Fj_names  = 'Fpv1 Fpv2 Fm'
    hj_names  = 'h h h'
    cj_names  = 'c2pv1 c2pv2 c2m'
    eta_i     = eta
    args      = 'c1m c1pv2 c1pv1  eta'
  [../]
  [./ACInterfacepv1]
    type = ACInterface
    variable = eta
    kappa_name = kappa
  [../]

  # Kernels for Allen-Cahn equation for eta_pv2
  [./deta_pv2_dt]
    type = TimeDerivative
    variable = eta
  [../]
  [./ACBulkFpv2]
    type = KKSMultiACBulkF
    variable  = eta
    Fj_names  = 'Fpv1 Fpv2 Fm'
    hj_names  = 'h h h'
    gi_name   = g
    eta_i     = eta
    wi        = 0.01
    args      = 'c1m c2m c1pv1 c2pv1 c1pv2 c2pv2 eta'
  [../]
  [./ACBulkCpv2_c1]
    type = KKSMultiACBulkC
    variable  = eta
    Fj_names  = 'Fpv1 Fpv2 Fm'
    hj_names  = 'h h h'
    cj_names  = 'c1pv1 c1pv2 c1m'
    eta_i     = eta
    args      = 'c2m c2pv2 c2pv1 eta'
  [../]
  [./ACBulkCpv2_c2]
    type = KKSMultiACBulkC
    variable  = eta
    Fj_names  = 'Fpv1 Fpv2 Fm'
    hj_names  = 'h h h'
    cj_names  = 'c2pv1 c2pv2 c2m'
    eta_i     = eta
    args      = 'c1m c1pv2 c1pv1 eta'
  [../]
  [./ACInterfacepv2]
    type = ACInterface
    variable = eta
    kappa_name = kappa
  [../]
[]

[AuxKernels]
  [./Energy_total]
    type = KKSMultiFreeEnergy
    Fj_names = 'Fpv1 Fpv2 Fm'
    hj_names = 'h h h'
    gj_names = 'g g g'
    variable = Energy
    w = 1
    interfacial_vars =  'eta  eta  eta'
    kappa_names =       'kappa kappa kappa'
  [../]
  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
    execute_on = timestep_end
    block = 0
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      ilu          nonzero'

  l_max_its = 50
  nl_max_its = 25
  l_tol = 1.0e-3
  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-8

  end_time = 14400

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

[Outputs]
  exodus = true
[]
