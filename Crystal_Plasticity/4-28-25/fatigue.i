[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [file]
     type = FileMeshGenerator
     file = Conc1_out.e-s202
     use_for_exodus_restart = true
   []
 []

 [Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]

  # order parameter 0
  [./eta0]
    initial_from_file_var = eta1
  [../]
  # order parameter 1
  [./eta1]
    initial_from_file_var = eta3
  [../]
[]



[BCs]
  [stressfree_boundary]
    #Applies the pressure
    type = Pressure
    boundary = top
    factor = 0.0 
    variable = disp_y
  []
  [dispy]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [dispx]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
 [right_pressure]
    type = Pressure
    variable = disp_x
    boundary = right
    function = cyclic_load
  []
[]

[Functions]
  [cyclic_load]
    type = ParsedFunction
    expression = -100*sin(2*pi*t/60) 
  []
[]

[UserObjects]
  [prop_read]
    type = PropertyReadFile
    prop_file_name = 'input_file.txt'
   # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 4
    read_type = voronoi
    nvoronoi = 1
    use_random_voronoi = true
    rand_seed = 67892
    rve_type = periodic
  []
[]


[AuxVariables]
  [bnds]
    order = FIRST
    family = LAGRANGE
  []
  [bounds_dummy]
    order = FIRST
    family = LAGRANGE
  []
  [accslip]
    order = CONSTANT
    family = MONOMIAL
  []
  [vonmises]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [fp_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [e_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [gss1]
    order = CONSTANT
    family = MONOMIAL
  []
[]


[Physics/SolidMechanics/QuasiStatic/all]
  strain = FINITE
  incremental = true
  add_variables = true
  strain_base_name = uncracked
[]

[AuxKernels]
  
  [vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    scalar_type = VonMisesStress
    execute_on = timestep_end
#   block = 0
  []
   [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  []
  [fp_xx]
    type = RankTwoAux
    variable = fp_xx
    rank_two_tensor = fp
    index_j = 0
    index_i = 0
    execute_on = timestep_end
#    block = 0
  []
  [e_xx]
    type = RankTwoAux
    variable = e_xx
    rank_two_tensor = lage
    index_j = 0
    index_i = 0
    execute_on = timestep_end
#    block = 0
  []
  [gss1]
    type = MaterialStdVectorAux
    variable = gss1
    property = gss
    index = 0
    execute_on = timestep_end
#    block = 0
  []
[]


[Materials]
  
  [elasticity_tensor_phase0]
    type = ComputeElasticityTensorCP
    C_ijkl = '2.906e5 1.87e5 1.87e5 2.906e5 1.87e5 2.906e5 1.142e5 1.142e5 1.142e5'
    fill_method = symmetric9
    base_name = phase0
  []
  [stress_phase0]
    type = ComputeMultipleCrystalPlasticityStress_abs
    crystal_plasticity_models = 'trial_xtalpl_phase0'
    tan_mod_type = exact
    rtol = 1e-08
    base_name = phase0
  []
  [trial_xtalpl_phase0]
    type = CrystalPlasticityKalidindiUpdate_abs
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
    crystal_lattice_type = FCC
    resistance_tol = 0.01
    r = 1.0             
    h = 6000            
    t_sat = 598.5        
    gss_a = 1.5         
    ao = 0.001           
    xm = 0.017             
    gss_initial = 600 
    base_name = phase0
  []
  [./strain_phase0]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y'
    base_name = phase0
    eigenstrain_names = eigenstrain1
  [../]
 [./eigen_strain0]
   type = ComputeEigenstrain
   base_name = phase0
   eigen_base = '1 1 0 0 0 0'
   prefactor = misfit
   eigenstrain_name = eigenstrain1
 [../]

  [elasticity_tensor_phase1]
    type = ComputeElasticityTensorCP
    C_ijkl = '2.721e5 1.69e5 1.69e5 2.721e5 1.69e5 2.721e5 1.31e5 1.31e5 1.31e5'
    fill_method = symmetric9
    base_name = phase1
  []
  [stress_phase1]
    type = ComputeMultipleCrystalPlasticityStress_abs
    crystal_plasticity_models = 'trial_xtalpl_phase1'
    tan_mod_type = exact
    rtol = 1e-08
    base_name = phase1
  []
  [trial_xtalpl_phase1]
    type = CrystalPlasticityKalidindiUpdate_abs
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
    crystal_lattice_type = FCC
    resistance_tol = 0.01
    r = 1.0             
    h = 6000            
    t_sat = 598.5        
    gss_a = 1.5         
    ao = 0.001           
    xm = 0.017             
    gss_initial = 465.5 
    base_name = phase1
  []
  [./strain_phase1]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y'
    base_name = phase1
  [../]

  
  # Switching functions for each phase
  [./h0]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta0
    all_etas = 'eta0 eta1'
    h_name = h0
  [../]
  [./h1]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta1
    all_etas = 'eta0 eta1'
    h_name = h1
  [../]
 
  # Generate the global stress from the phase stresses
  [./global_stress]
    type = MultiPhaseStressMaterial
    phase_base = 'phase0 phase1'
    h          = 'h0     h1'
  [../]

  [elasticity_tensor]
    type = ComputeElasticityTensorCP
    C_ijkl = '2.04e5 1.37e5 1.37e5 2.044e5 1.37e5 2.04e5 1.26e5 1.26e5 1.26e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
    output_properties = 'Euler_angles'
    outputs = exodus
    base_name = uncracked
  []
 
[]

[Kernels]
  [./eta0_dt]
    type = TimeDerivative
    variable = eta0
  [../]
  [./eta1_dt]
    type = TimeDerivative
    variable = eta1
  [../]
[]



[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu            superlu_dist          vinewtonrsls'  

#  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
#  petsc_options_value = ' asm      2              lu            gmres     200'
  l_max_its = 20
  nl_max_its = 10
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-9

  end_time = 50

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 5e-4
    cutback_factor = 0.75
    growth_factor = 1.2
    optimal_iterations = 20
  [../]

[]

# Precondition using handcoded off-diagonal terms

[Preconditioning]
  [full]
    type = SMP
    full = true
  []
[]


[Postprocessors]
[accslip]
    type = ElementAverageValue
    variable = accslip
    block = 'ANY_BLOCK_ID 0'
    []
 [vonmises]
    type = ElementAverageValue
    variable = vonmises
    block = 'ANY_BLOCK_ID 0'
  []
  [stress_xx]
    type = ElementAverageValue
    variable = stress_xx
    block = 'ANY_BLOCK_ID 0'
  []
  [fp_xx]
    type = ElementAverageValue
    variable = fp_xx
    block = 'ANY_BLOCK_ID 0'
  []
  [e_xx]
    type = ElementAverageValue
    variable = e_xx
    block = 'ANY_BLOCK_ID 0'
  []
  [gss1]
    type = ElementAverageValue
    variable = gss1
    block = 'ANY_BLOCK_ID 0'
  []
[]

[Outputs]
  exodus = true
  time_step_interval = 2
  [console]
     type = Console
     time_step_interval = 2
  []
  [csv]
     type = CSV
     time_step_interval = 2
  []
[]