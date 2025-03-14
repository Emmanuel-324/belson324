#
# Fig. 6 input for 10.1016/j.commatsci.2017.02.017
# D. Schwen et al./Computational Materials Science 132 (2017) 36-45
# Three phase interface simulation demonstrating the interfacial stability
# w.r.t. formation of a tspurious third phase
#

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 120
  ny = 120
  nz = 0
  xmin = 0
  xmax = 40
  ymin = 0
  ymax = 40
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[Variables]
   # concentration
   [./c]
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

  # phase concentration 1
  [./c1]
    order = FIRST
    family = LAGRANGE
#    initial_condition = 0.3333
  [../]

  # phase concentration 2
  [./c2]
    order = FIRST
    family = LAGRANGE
#    initial_condition = 0.3333
  [../]

  # phase concentration 3
  [./c3]
    order = FIRST
    family = LAGRANGE
#    initial_condition = 0.2000
  [../]
  # Lagrange multiplier
  [./lambda]
    initial_condition = 0.0
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

[AuxVariables]
  [./T]
    [./InitialCondition]
      type = FunctionIC
      function = 'x-10'
    [../]
  [../]
[]


[Functions]
  [./ic_func_eta1]
    type = ParsedFunction
    expression = '0.5*(1.0+tanh((x-10)/sqrt(2.0))) * 0.5*(1.0+tanh((y-10)/sqrt(2.0)))'
  [../]
  [./ic_func_eta2]
    type = ParsedFunction
    expression = '0.5*(1.0-tanh((x-10)/sqrt(2.0)))'
  [../]
  [./ic_func_eta3]
    type = ParsedFunction
    expression = '1 - 0.5*(1.0-tanh((x-10)/sqrt(2.0)))
              - 0.5*(1.0+tanh((x-10)/sqrt(2.0))) * 0.5*(1.0+tanh((y-10)/sqrt(2.0)))'
  [../]
  [./ic_func_c]
    type = ParsedFunction
    expression = '0.5 * 0.5*(1.0-tanh((x-10)/sqrt(2.0)))
              + 0.4 * 0.5*(1.0+tanh((x-10)/sqrt(2.0))) * 0.5*(1.0+tanh((y-10)/sqrt(2.0)))
              + 0.8 * (1 - 0.5*(1.0-tanh((x-10)/sqrt(2.0)))
                        - 0.5*(1.0+tanh((x-10)/sqrt(2.0))) * 0.5*(1.0+tanh((y-10)/sqrt(2.0))))'
  [../]
[]

[ICs]
  [./eta1]
    variable = eta1
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 192
  [../]
  [./eta2]
    variable = eta2
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 389	
  [../]
  [./c]
    variable = c
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
    property_name = F1
    coupled_variables = 'c1'
    expression = '100.0*(c1-0.3111)^2'
  [../]
  [./f2]
    type = DerivativeParsedMaterial
    property_name = F2
    coupled_variables = 'c2'
     expression = '100.0*(c2-0.34)^2'
  [../]
  [./f3]
    type = DerivativeParsedMaterial
    property_name = F3
    coupled_variables = 'c3'
    expression = '5.0*(c3-0.4)^2'
  [../]

  # Switching functions for each phase
  # h1(eta1, eta2, eta3)
  [./h1]
    type = SwitchingFunction3PhaseMaterial
    eta_i = eta1
    eta_j = eta2
    eta_k = eta3
    f_name = h1
  [../]
  # h2(eta1, eta2, eta3)
  [./h2]
    type = SwitchingFunction3PhaseMaterial
    eta_i = eta2
    eta_j = eta3
    eta_k = eta1
    f_name = h2
  [../]
  # h3(eta1, eta2, eta3)
  [./h3]
    type = SwitchingFunction3PhaseMaterial
    eta_i = eta3
    eta_j = eta1
    eta_k = eta2
    f_name = h3
  [../]

  # Coefficients for diffusion equation
  [./Dh1]
    type = DerivativeParsedMaterial
    material_property_names = 'D h1'
    expression = D*h1
    property_name = Dh1
  [../]
  [./Dh2]
    type = DerivativeParsedMaterial
    material_property_names = 'D h2'
    expression = D*h2
    property_name = Dh2
  [../]
  [./Dh3]
    type = DerivativeParsedMaterial
    material_property_names = 'D h3'
    expression = D*h3
    property_name = Dh3
  [../]

  # Barrier functions for each phase
  [./g1]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta1
    function_name = g1
  [../]
  [./g2]
    type = BarrierFunctionMaterial
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

[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
  #Kernels for diffusion equation
  [./diff_time]
    type = TimeDerivative
    variable = c
  [../]
  [./diff_c1]
    type = MatDiffusion
    variable = c
    diffusivity = Dh1
    v = c1
  [../]
  [./diff_c2]
    type = MatDiffusion
    variable = c
    diffusivity = Dh2
    v = c2
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
    args      = 'c1 c2 c3 eta2 eta3'
  [../]
  [./ACBulkC1]
    type = KKSMultiACBulkC
    variable  = eta1
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta1
    args      = 'eta2 eta3'
  [../]
  [./ACInterface1]
    type = ACInterface
    variable = eta1
    kappa_name = kappa
  [../]
  [./multipler1]
    type = MatReaction
    variable = eta1
    v = lambda
    reaction_rate = L
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
    args      = 'c1 c2 c3 eta1 eta3'
  [../]
  [./ACBulkC2]
    type = KKSMultiACBulkC
    variable  = eta2
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta2
    args      = 'eta1 eta3'
  [../]
  [./ACInterface2]
    type = ACInterface
    variable = eta2
    kappa_name = kappa
  [../]
  [./multipler2]
    type = MatReaction
    variable = eta2
    v = lambda
    reaction_rate = L
  [../]

  # Kernels for the Lagrange multiplier equation
  [./mult_lambda]
    type = MatReaction
    variable = lambda
    reaction_rate = 3
  [../]
  [./mult_ACBulkF_1]
    type = KKSMultiACBulkF
    variable  = lambda
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g1
    eta_i     = eta1
    wi        = 0.01
    mob_name  = 1
    args      = 'c1 c2 c3 eta2 eta3'
  [../]
  [./mult_ACBulkC_1]
    type = KKSMultiACBulkC
    variable  = lambda
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta1
    args      = 'eta2 eta3'
    mob_name  = 1
  [../]
  [./mult_CoupledACint_1]
    type = SimpleCoupledACInterface
    variable = lambda
    v = eta1
    kappa_name = kappa
    mob_name = 1
  [../]
  [./mult_ACBulkF_2]
    type = KKSMultiACBulkF
    variable  = lambda
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g2
    eta_i     = eta2
    wi        = 0.01
    mob_name  = 1
    args      = 'c1 c2 c3 eta1 eta3'
  [../]
  [./mult_ACBulkC_2]
    type = KKSMultiACBulkC
    variable  = lambda
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta2
    args      = 'eta1 eta3'
    mob_name  = 1
  [../]
  [./mult_CoupledACint_2]
    type = SimpleCoupledACInterface
    variable = lambda
    v = eta2
    kappa_name = kappa
    mob_name = 1
  [../]
  [./mult_ACBulkF_3]
    type = KKSMultiACBulkF
    variable  = lambda
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    gi_name   = g3
    eta_i     = eta3
    wi        = 0.01
    mob_name  = 1
    args      = 'c1 c2 c3 eta1 eta2'
  [../]
  [./mult_ACBulkC_3]
    type = KKSMultiACBulkC
    variable  = lambda
    Fj_names  = 'F1 F2 F3'
    hj_names  = 'h1 h2 h3'
    cj_names  = 'c1 c2 c3'
    eta_i     = eta3
    args      = 'eta1 eta2'
    mob_name  = 1
  [../]
  [./mult_CoupledACint_3]
    type = SimpleCoupledACInterface
    variable = lambda
    v = eta3
    kappa_name = kappa
    mob_name = 1
  [../]

  # Kernels for constraint equation eta1 + eta2 + eta3 = 1
  # eta3 is the nonlinear variable for the constraint equation
  [./eta3reaction]
    type = MatReaction
    variable = eta3
    reaction_rate = 1
  [../]
  [./eta1reaction]
    type = MatReaction
    variable = eta3
    v = eta1
    reaction_rate = 1
  [../]
  [./eta2reaction]
    type = MatReaction
    variable = eta3
    v = eta2
    reaction_rate = 1
  [../]
  [./one]
    type = BodyForce
    variable = eta3
    value = -1.0
  [../]

  # Phase concentration constraints
  [./chempot12]
    type = KKSPhaseChemicalPotential
    variable = c1
    cb       = c2
    fa_name  = F1
    fb_name  = F2
  [../]
  [./chempot23]
    type = KKSPhaseChemicalPotential
    variable = c2
    cb       = c3
    fa_name  = F2
    fb_name  = F3
  [../]
  [./phaseconcentration]
    type = KKSMultiPhaseConcentration
    variable = c3
    cj = 'c1 c2 c3'
    hj_names = 'h1 h2 h3'
    etas = 'eta1 eta2 eta3'
    c = c
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
  [./gr1area]
      type = ElementIntegralVariablePostprocessor_new2
      variable = eta1
      execute_on = 'initial timestep_end'
  [../]
    [./gr2area]
      type = ElementIntegralVariablePostprocessor_new2
      variable = eta2
      execute_on = 'initial timestep_end'
  [../]
    [./gr3area]
      type = ElementIntegralVariablePostprocessor_new2
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

#[VectorPostprocessors]
#  [./c]
#    type =  LineValueSampler
#    start_point = '-25 0 0'
#    end_point = '25 0 0'
#    variable = c
#    num_points = 151
#    sort_by =  id
#    execute_on = timestep_end
#  [../]
#  [./eta1]
#    type =  LineValueSampler
#    start_point = '-25 0 0'
#    end_point = '25 0 0'
#    variable = eta1
#    num_points = 151
#    sort_by =  id
#    execute_on = timestep_end
#  [../]
#  [./eta2]
#    type =  LineValueSampler
#    start_point = '-25 0 0'
#    end_point = '25 0 0'
#    variable = eta2
#    num_points = 151
#    sort_by =  id
#    execute_on = timestep_end
#  [../]
#  [./eta3]
#    type =  LineValueSampler
#    start_point = '-25 0 0'
#    end_point = '25 0 0'
#    variable = eta3
#    num_points = 151
#    sort_by =  id
#    execute_on = timestep_end
#  [../]
#[]