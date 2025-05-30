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
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = eta1
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_lower_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = eta1
    bound_type = lower
    bound_value = -1
  [../]
  [./eta2_upper_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = eta2
    bound_type = upper
    bound_value = 1
  [../]
  [./eta2_lower_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = eta2
    bound_type = lower
    bound_value = -1
  [../]
[]


[Variables]
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

# Niobium (Ni) solute concentration
[./c_Ni]
  order = FIRST
  family = LAGRANGE
[../]
  
  # order parameter 1
  [./eta1]
    order = FIRST
    family = LAGRANGE
  [../]

  # order parameter 2
  [./eta2]
    order = FIRST
    family = LAGRANGE
  [../]

  # order parameter 3
  [./eta3]
    order = FIRST
    family = LAGRANGE
#    initial_condition = 0.0
  [../]

  # Gamma (γ) phase solute concentrations
  [./c_Al_gamma]
    order = FIRST
    family = LAGRANGE
   initial_condition = 0.1
  [../]

 
  [./c_Nb_gamma]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.05
  [../]

  # Gamma Prime (γ') phase solute concentrations
  [./c_Al_gammaP]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.8
  [../]

  # Solid phase solute 2 concentration
  [./c_Nb_gammaP]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.1
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
  [./ic_func_c]
    type = ParsedFunction
    value = 0.5+0.01*(cos(1.05*x)*cos(1.1*y)+(cos(1.3*x)*cos(0.87*y))^2+cos(0.25*x-1.5*y)*cos(0.7*x-0.2*y))
  [../]
  [./bc_func]
    type = ParsedFunction
    value = sin(alpha*pi*x)
    vars = alpha
    vals = 16
  [../]
  [./disp_func]
    type = ParsedFunction
    value = 'if(t<50,6e-3*t,0.3)'
  [../]
  [./press_func]
    type = ParsedFunction
    value = '1'
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
    seed = 89	
  [../]
  [./c_Nb]
    variable = c_Nb
    type = RandomIC
    min = 0.24375	
    max = 0.25625
    seed = 89	
  [../]
[]


[Materials]
  # simple toy free energies
  [./f1]
    type = DerivativeParsedMaterial
    f_name = fc_1
    args = 'c_Al_gamma c_Nb_gamma'
    function = '1e5 * (4e5 * (c_Al_gamma - 0.187)^2 + 9.5e5 * (c_Nb_gamma- 0.0157)^2)'
  [../]
  # Elastic energy of the phase 1
  [./elastic_free_energy_1]
    type = ElasticEnergyMaterial
    base_name = phase1
    f_name = fe_1
    args = ' '
  [../]
  # Total free energy of the phase 1
  [./Total_energy_1]
    type = DerivativeSumMaterial
    f_name = F1
    sum_materials = 'fc_1 fe_1'
    args = 'c_Al_gamma c_Nb_gamma'
  [../]

  [./f2]
    type = DerivativeParsedMaterial
    f_name = fc_2
    args = 'c_Al_gammaP c_Nb_gammaP'
    function = '1e5 * (4e5 * (c_Al_gammaP - 0.000727)^2 + 9.5e5 * (c_Nb_gammaP - 0.196)^2) + 1.5485e8'
  [../]
  # Elastic energy of the phase 2
  [./elastic_free_energy_2]
    type = ElasticEnergyMaterial
    base_name = phase2
    f_name = fe_2
    args = ' '
  [../]
  # Total free energy of the phase 2
  [./Total_energy_2]
    type = DerivativeSumMaterial
    f_name = F2
    sum_materials = 'fc_2 fe_2'
    args ='c_Al_gammaP c_Nb_gammaP'
  [../]

  [./f3]
    type = DerivativeParsedMaterial
    f_name = fc_3
    args = 'c3'
    function = '5.0*(c3-0.20)^2'
  [../]
  # Elastic energy of the phase 3
  [./elastic_free_energy_3]
    type = ElasticEnergyMaterial
    base_name = phase3
    f_name = fe_3
    args = ' '
  [../]
  # Total free energy of the phase 3
  [./Total_energy_3]
    type = DerivativeSumMaterial
    f_name = F3
    sum_materials = 'fc_3 fe_3'
    args = 'c3'
  [../]

  # Switching functions for each phase
  # h1(eta1, eta2, eta3)
  [./h1]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta1
    all_etas = 'eta1 eta2 eta3'
    h_name = h1
  [../]
  # h2(eta1, eta2, eta3)
  [./h2]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta2
    all_etas = 'eta1 eta2 eta3'
    h_name = h2
  [../]
  # h3(eta1, eta2, eta3)
  [./h3]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta3
    all_etas = 'eta1 eta2 eta3'
    h_name = h3
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
  [./Dh3]
    type = DerivativeParsedMaterial
    material_property_names = 'D h3'
    function = D*h3
    f_name = Dh3
  [../]

  # Barrier functions for each phase
  [./g1]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta1
    function_name = g1
  [../]
  [./g2]
    type = BarrierFunctionMaterial_abs
    g_order = SIMPLE
    eta = eta2
    function_name = g2
  [../]
  [./g3]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta3
    function_name = g3
  [../]

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'L    kappa  D  misfit  W'
    prop_values = '0.3  0.01   1  0.005  0.01'
  [../]

  #Mechanical properties
  [./Stiffness_phase1]
    type = ComputeElasticityTensor
    C_ijkl = '1982 534 496 1606 605 1788 985 1056 388'    
    base_name = phase1
    fill_method = symmetric9
  [../]
  [./Stiffness_phase2]
    type = ComputeElasticityTensor
    C_ijkl = '1606 534 605 1982 496 1788 1056 985 388'
    base_name = phase2
    fill_method = symmetric9
  [../]
  [./Stiffness_phase3]
    type = ComputeElasticityTensor
    C_ijkl = '2008.5 448.9 918.5 2008.5 918.5 1538.7 779.8 779.8 309.9'    
    base_name = phase3
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
  [./stress_phase3]
    type = ComputeLinearElasticStress
    base_name = phase3
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
  [./strain_phase3]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phase3
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
    phase_base = 'phase1 phase2 phase3'
    h          = 'h1     h2     h3'
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

  #Kernels for diffusion equation
  [./diff_time]
    type = TimeDerivative
    variable = c_Al
  [../]
  [./diff_time1]
    type = TimeDerivative
    variable = c_Nb
  [../]
  [./diff_c_Al]
    type = MatDiffusion
    variable = c_Al
    diffusivity = Dh1
    v = c_Al
  [../]
  [./diff_c_Nb]
    type = MatDiffusion
    variable = c_Nb
    diffusivity = Dh2
    v = c_Nb
  [../]
  [./diff_c3]
    type = MatDiffusion
    variable = c
    diffusivity = Dh3
    v = c3
  [../]

  # Kernels for Allen-Cahn equation for eta1
  [./deta1dt]
    type = TimeDerivative
    variable = eta1
  [../]
  [./ACBulkF1]
    type = KKSMultiACBulkF
    variable  = eta1
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g1
    eta_i     = eta1
    wi        = 0.01
    args      = 'c_Al c_Nb c3 eta2 eta3'
  [../]
  [./ACBulkC1]
    type = KKSMultiACBulkC
    variable  = eta1
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c_Al c_Nb c3'
    eta_i     = eta1
    args      = 'eta2 eta3'
  [../]
  [./ACInterface1]
    type = ACInterface
    variable = eta1
    kappa_name = kappa
  [../]

  # Kernels for Allen-Cahn equation for eta2
  [./deta2dt]
    type = TimeDerivative
    variable = eta2
  [../]
  [./ACBulkF2]
    type = KKSMultiACBulkF
    variable  = eta2
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g2
    eta_i     = eta2
    wi        = 0.01
    args      = 'c_Al c_Nb c3 eta1 eta3'
  [../]
  [./ACBulkC2]
    type = KKSMultiACBulkC
    variable  = eta2
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c_Al c_Nb c3'
    eta_i     = eta2
    args      = 'eta1 eta3'
  [../]
  [./ACInterface2]
    type = ACInterface
    variable = eta2
    kappa_name = kappa
  [../]

  # Kernels for constraint equation eta1 + eta2 + eta3 = 1
  # eta3 is the nonlinear variable for the constraint equation
  [./eta3reaction]
    type = MatReaction
    variable = eta3
    mob_name = 1
  [../]
  [./eta1reaction]
    type = MatReaction_abscouple
    variable = eta3
    v = eta1
    mob_name = 1
  [../]
  [./eta2reaction]
    type = MatReaction_abscouple
    variable = eta3
    v = eta2
    mob_name = 1
  [../]
  [./one]
    type = BodyForce
    variable = eta3
    value = -1.0
  [../]


  # Phase concentration constraints
  [./chempot12]
    type = KKSPhaseChemicalPotential
    variable = c_Al
    cb       = c_Nb
    fa_name  = F1
    fb_name  = F2
  [../]
  [./chempot23]
    type = KKSPhaseChemicalPotential
    variable = c_Nb
    cb       = c3
    fa_name  = F2
    fb_name  = F3
  [../]
  [./phaseconcentration]
    type = KKSMultiPhaseConcentration
    variable = c3
    cj = 'c_Al c_Nb c3'
    hj_names = 'h1 h2 h3'
    etas = 'eta1 eta2 eta3'
    c = c
  [../]
[]

[AuxKernels]
  [./Energy_total]
    type = KKSMultiFreeEnergy
    Fj_names = 'F1 F2 F3'
    hj_names = 'h1 h2 h3'
    gj_names = 'g1 g2 g3'
    variable = Energy
    w = 1
    interfacial_vars =  'eta1  eta2  eta3'
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
    [./gr3area]
      type = ElementIntegralVariablePostprocessor
      variable = eta3
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