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
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Bounds]
  [./eta_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta1
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta1
    bound_type = lower
    bound_value = -1
  [../]
  [./eta2_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta2
    bound_type = upper
    bound_value = 1
  [../]
  [./eta2_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta2
    bound_type = lower
    bound_value = -1
  [../]
[]

[BCs]
  [./u_right_pull]
    type = Pressure
    displacements = 'disp_x disp_y'
    variable = disp_x
    boundary = right
    factor = 0
  [../]
  [./c_Al]
    type =  NeumannBC
    variable = 'c_Al'
    boundary = 'left right top bottom'
    value = 0
  [../]	
  [./c_Nb]
    type =  NeumannBC
    variable = 'c_Nb'
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

[Variables]
  # order parameter
  [./eta1]
    order = FIRST
    family = LAGRANGE
  [../]
  
  # order parameter
  [./eta2]
    order = FIRST
    family = LAGRANGE
  [../]

# Aluminum (Al) solute concentration  
  [./c_Al]
    order = FIRST
    family = LAGRANGE
    
  [../]

 # Niobium (Nb) solute concentration
  [./c_Nb]
    order = FIRST
    family = LAGRANGE
    
  [../]

  # chemical potential solute 1
  [./w_Al]
    order = FIRST
    family = LAGRANGE
  [../]

  # chemical potential solute 2
  [./w_Nb]
    order = FIRST
    family = LAGRANGE
  [../]

# Gamma (γ) phase solute concentrations
  [./c_Al_gamma]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0161
  [../]

 
  [./c_Nb_gamma]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.00723
  [../]

  # Gamma Prime (γ') phase solute concentrations
  [./c_Al_gammaP]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.187
  [../]

  # Solid phase solute 2 concentration
  [./c_Nb_gammaP]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0157
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
  [./ic_func_eta1]
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
  [./bc_func]
    type = ParsedFunction
    expression = sin(alpha*pi*x)
    vars = alpha
    vals = 16
  [../]
[]

[ICs]
  [./eta1]
    variable = eta1
    type = RandomIC
    min = -0.1625
    max = 0.1625
    seed = 192
  [../]
  [./eta2]
    variable = eta2
    type = RandomIC
    min = -0.1625
    max = 0.1625
    seed = 389	
  [../]
  [./c_Al]
    variable = c_Al
    type = RandomIC
    min = 0.24375
    max = 0.25625
    seed = 389	
  [../]
  [./c_Nb]
    variable = c_Nb
    type = RandomIC
    min = 0.24375
    max = 0.25625
    seed = 389	
  [../]
[]

[Materials]
  # Free energy of the gamma
  [./f_gamma]
    type = DerivativeParsedMaterial
    property_name = f_gamma
    coupled_variables = 'c_Al_gamma c_Nb_gamma'
    expression = '(0.0161-c_Al_gamma)^2+(0.00723-c_Nb_gamma)^2'
  [../]
   # Elastic energy of the gamma
   [./elastic_free_energy_1]
    type = ElasticEnergyMaterial
    base_name = phase1
    f_name = fe_gamma
    coupled_variables = ' '
  [../]
  # Total free energy of the gamma
  [./Total_energy_1]
    type = DerivativeSumMaterial
    f_name = F1
    sum_materials = 'f_gamma fe_gamma'
    coupled_variables = 'c_Al_gamma c_Nb_gamma'
  [../]
  # Free energy of the gammaP
  [./f_gammaP]
    type = DerivativeParsedMaterial
    property_name = f_gammaP
    coupled_variables = 'c_Al_gammaP c_Nb_gammaP'
    expression =  '(0.187-c_Al_gammaP)^2+(0.0157-c_Nb_gammaP)^2'
  [../]
   # Elastic energy of the phase 2
   [./elastic_free_energy_2]
    type = ElasticEnergyMaterial
    base_name = phase2
    f_name = fe_gammaP
    coupled_variables = ' '
  [../]
   # Total free energy of the gammaP
  [./Total_energy_2]
    type = DerivativeSumMaterial
    f_name = F2
    sum_materials = 'f_gammaP fe_gammaP'
    coupled_variables = 'c_Al_gammaP c_Nb_gammaP'
  [../]
  # h(eta1)
  [./h1]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta1
    all_etas = 'eta1 eta2'
    h_name = h1
  [../]
  [./h2]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta2
    all_etas = 'eta1 eta2 '
    h_name = h2
  [../]
# Coefficients for diffusion equation
[./Dh1]
  type = DerivativeParsedMaterial
  material_property_names = 'D h1'
  function = D*h1
  f_name = Dh1
[../]
[./Dh2]
  type = DerivativeParsedMaterial
  material_property_names = 'D h2'
  function = D*h2
  f_name = Dh2
[../]
  # g(eta1)
  [./g_eta1]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta1
    function_name = g1
  [../]
  # g(eta2)
  [./g_eta2]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta2
    function_name = g2
  [../]

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'M       L     kappa  eps_sq    misfit D'
    prop_values = '1e-3   0.3     5e-9     1.0     0.005  1'
  [../]
   #Mechanical properties
   [./Stiffness_phase1]
    type = ComputeElasticityTensor
    C_ijkl = '2.721e5 1.69e5 1.69e5 2.721e5 1.69e5 2.721e5 1.31e5 1.31e5 1.31e5'    
    base_name = phase1
    fill_method = symmetric9
  [../]
  [./Stiffness_phase2]
    type = ComputeElasticityTensor
    C_ijkl = '2.721e5 1.69e5 1.69e5 2.721e5 1.69e5 2.721e5 1.31e5 1.31e5 1.31e5'
    base_name = phase2
    fill_method = symmetric9
  [../]
 
  [./stress_phase1]
    type = ComputeLinearElasticStress
    base_name = phase1
  [../]
  [./stress_phase2]
    type = ComputeLinearElasticStress
    base_name = phase2
  [../]
 

  [./strain_phase1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phase1
    eigenstrain_names = eigenstrain1
  [../]
  [./strain_phase2]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phase2
    eigenstrain_names = eigenstrain2
  [../]
 
  [./eigen_strain1]
    type = ComputeEigenstrain
    base_name = phase1
    eigen_base = '1 0 0 0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrain1
  [../]

  [./eigen_strain2]
    type = ComputeEigenstrain
    base_name = phase2
    eigen_base = '0 1 0 0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrain2
  [../]

  # Generate the global stress from the phase stresses
  [./global_stress]
    type = MultiPhaseStressMaterial
    phase_base = 'phase1 phase2'
    h          = 'h1     h2 '
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
  [./diff_time_Al]
    type = TimeDerivative
    variable = c_Al
  [../]
  [./diff_time_Nb]
    type = TimeDerivative
    variable = c_Nb
  [../]
  [./diff_c_Al_gamma]
    type = MatDiffusion
    variable = c_Al
    diffusivity = Dh1
    v = c_Al_gamma
  [../]
  [./diff_c_Nb_gamma]
    type = MatDiffusion
    variable = c_Nb
    diffusivity = Dh1
    v = c_Nb_gamma
  [../]
  [./diff_c_Al_gammaP]
    type = MatDiffusion
    variable = c_Al
    diffusivity = Dh2
    v = c_Al_gammaP
  [../]

  [./diff_c_Nb_gammaP]
    type = MatDiffusion
    variable = c_Nb
    diffusivity = Dh2
    v = c_Nb_gammaP
  [../]

  # enforce c_Al = (1-h(eta1))*c_Al_gamma + h(eta1)*c_Al_gammaP
  [./PhaseConc_c_Al_phase1]
    type = KKSPhaseConcentration
    ca       = c_Al_gamma
    variable = c_Al_gammaP
    c        = c_Al
    eta      = eta1
    h_name = h1
  [../]
  
  # enforce c2 = (1-h(eta1))*c_Nb_gamma + h(eta1)*c_Nb_gammaP
  [./PhaseConc_c_Nb_phase1]
    type = KKSPhaseConcentration
    ca       = c_Nb_gamma
    variable = c_Nb_gammaP
    c        = c_Nb
    eta      = eta1
    h_name = h1
  [../]

 # enforce c1 = (1-h(eta2))*c_Al_gamma + h(eta2)*c_Al_gammaP
 [./PhaseConc_c_Al_phase2]
  type = KKSPhaseConcentration
  ca       = c_Al_gamma
  variable = c_Al_gammaP
  c        = c_Al
  eta      = eta2
  h_name = h2
[../]

  # enforce c2 = (1-h(eta2))*c_Nb_gamma + h(eta2)*c_Nb_gammaP
  [./PhaseConc_c_Nb_phase2]
    type = KKSPhaseConcentration
    ca       = c_Nb_gamma
    variable = c_Nb_gammaP
    c        = c_Nb
    eta      = eta2
    h_name = h2
  [../]

  # enforce pointwise equality of chemical potentials
  [./ChemPotSolute1]
    type = KKSPhaseChemicalPotential
    variable = c_Al_gamma
    cb       = c_Al_gammaP
    fa_name  = f_gamma
    fb_name  = f_gammaP
    args_a   = 'c_Nb_gamma'
    args_b   = 'c_Nb_gammaP'
  [../]
  [./ChemPotSolute2]
    type = KKSPhaseChemicalPotential
    variable = c_Nb_gamma
    cb       = c_Nb_gammaP
    fa_name  = f_gamma
    fb_name  = f_gammaP
    args_a   = 'c_Al_gamma'
    args_b   = 'c_Al_gammaP'

  [../]

  #
  # Cahn-Hilliard Equations
  #
  [./CHBulk1]
    type = KKSSplitCHCRes
    variable = c_Al
    ca       = c_Al_gamma
    fa_name  = f_gamma
    w        = w_Al
    args_a   = 'c_Nb_gamma'
  [../]
  [./CHBulk2]
    type = KKSSplitCHCRes
    variable = c_Nb
    ca       = c_Nb_gamma
    fa_name  = f_gamma
    w        = w_Nb
    args_a   = 'c_Al_gamma'
  [../]

  [./dc_Aldt]
    type = CoupledTimeDerivative
    variable = w_Al
    v = c_Al
  [../]
  [./dc_Nbdt]
    type = CoupledTimeDerivative
    variable = w_Nb
    v = c_Nb
  [../]

  [./w_Alkernel]
    type = SplitCHWRes
    mob_name = M
    variable = w_Al
  [../]
  [./w_Nbkernel]
    type = SplitCHWRes
    mob_name = M
    variable = w_Nb
  [../]

  #
  # Allen-Cahn Equation
  #
  [./ACBulkF1]
    type = KKSACBulkF
    variable = eta1
    fa_name  = f_gamma
    fb_name  = f_gammaP
    h_name = h1
    g_name = g1
    w        = 1e5
    coupled_variables = 'c_Al_gamma c_Al_gammaP c_Nb_gamma c_Nb_gammaP'
  [../]
  [./ACBulkF2]
    type = KKSACBulkF
    variable = eta2
    fa_name  = f_gamma
    fb_name  = f_gammaP
    h_name = h2
    g_name = g2
    w        = 1e5
    coupled_variables = 'c_Al_gamma c_Al_gammaP c_Nb_gamma c_Nb_gammaP'
  [../]
  [./ACBulkC1_eta1]
    type = KKSACBulkC
    variable = eta1
    ca       = c_Al_gamma
    cb       = c_Al_gammaP
    h_name = h1
    fa_name  = f_gamma
    coupled_variables     = 'c_Nb_gamma'
  [../]
  [./ACBulkC1_eta2]
    type = KKSACBulkC
    variable = eta2
    ca       = c_Al_gamma
    cb       = c_Al_gammaP
    h_name = h2
    fa_name  = f_gamma
    coupled_variables     = 'c_Nb_gamma'
  [../]
  [./ACBulkC2_eta1]
    type = KKSACBulkC
    variable = eta1
    ca       = c_Nb_gamma
    cb       = c_Nb_gammaP
    h_name = h1
    fa_name  = f_gamma
    coupled_variables     = 'c_Al_gamma'
  [../]
  [./ACBulkC2_eta2]
    type = KKSACBulkC
    variable = eta2
    ca       = c_Nb_gamma
    cb       = c_Nb_gammaP
    h_name = h2
    fa_name  = f_gamma
    coupled_variables     = 'c_Al_gamma'
  [../]
  [./ACInterface1]
    type = ACInterface
    variable = eta1
    kappa_name = kappa
  [../]
  [./ACInterface2]
    type = ACInterface
    variable = eta2
    kappa_name = kappa
  [../]
  [./deta1dt]
    type = TimeDerivative
    variable = eta1
  [../]
  [./deta2dt]
      type = TimeDerivative
      variable = eta2
  [../]
[]

[AuxKernels]
  [./Energy_total]
    type = KKSMultiFreeEnergy
    Fj_names = 'F1 F2'
    hj_names = 'h1 h2'
    gj_names = 'g1 g2'
    variable = Energy
    w = 1
    interfacial_vars =  'eta1  eta2 '
    kappa_names =       'kappa kappa'
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
#  [./e_xx]
#    type = RankTwoAux
#    variable = e_xx
#    rank_two_tensor = phase1_lage
#    index_j = 0
#    index_i = 0
#    execute_on = timestep_end
#    block = 0
#  [../]
[]
[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu            mumps            vinewtonrsls'
  line_search = none
  l_max_its = 50
  nl_max_its = 25
  l_tol = 1.0e-3
  nl_rel_tol = 1.0e-7
  nl_abs_tol = 1.0e-9

  end_time = 100

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 5e-6
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
    variable = eta1
    function = bc_func
  [../]
   [./gr1area]
     type = ElementIntegralVariablePostprocessor
     variable = eta1
     execute_on = 'initial timestep_end'
 [../]
   [./gr2area]
     type = ElementIntegralVariablePostprocessor
     variable = eta2
     execute_on = 'initial timestep_end'
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
     interval = 10
  [../]
[]