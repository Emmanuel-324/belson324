[GlobalParams]
  order = FIRST
  family = LAGRANGE
  displacements = 'disp_x disp_y'
  temperature = temperature
  out_of_plane_strain = strain_zz
[]

[Mesh]
  use_displaced_mesh = true
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
  [./strain_zz]
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
[AuxVariables]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [temperature]
    order = FIRST
    family = LAGRANGE
  []
  [./nl_strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [eth_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [eth_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [eth_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  
  [fth_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [fth_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [fth_zz]
    order = CONSTANT
    family = MONOMIAL
  []
[]
[Physics/SolidMechanics/QuasiStatic]
  [plane_stress]
    planar_formulation = WEAK_PLANE_STRESS
    strain = FINITE
    add_variables = true
    generate_output = 'stress_xx stress_xy stress_yy stress_zz strain_xx strain_xy strain_yy'
    eigenstrain_names = deformation_gradient
  []
[]

[AuxKernels]
  [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    scalar_type = VonMisesStress
    execute_on = timestep_end
  [../]
  [temperature]
    type = FunctionAux
    variable = temperature
    function = '298'
    execute_on = timestep_begin
  []
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = nl_strain_zz
    index_i = 2
    index_j = 2
  [../]
  [eth_xx]
    type = RankTwoAux
    variable = eth_xx
    rank_two_tensor = phase0_thermal_eigenstrain
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  []
  [eth_yy]
    type = RankTwoAux
    variable = eth_yy
    rank_two_tensor = phase0_thermal_eigenstrain
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  []
  [eth_zz]
    type = RankTwoAux
    variable = eth_zz
    rank_two_tensor = phase0_thermal_eigenstrain
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []
  [fth_xx]
    type = RankTwoAux
    variable = fth_xx
    rank_two_tensor = phase0_thermal_deformation_gradient
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  []
  [fth_yy]
    type = RankTwoAux
    variable = fth_yy
    rank_two_tensor = phase0_thermal_deformation_gradient
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  []
  [fth_zz]
    type = RankTwoAux
    variable = fth_zz
    rank_two_tensor = phase0_thermal_deformation_gradient
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []
  
[]

[BCs]
  [symmy]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [symmx]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  
  [tdisp]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = '0.1*t'
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
    eigenstrain_names = 'thermal_eigenstrain'
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
  [thermal_eigenstrain]
    type = ComputeCrystalPlasticityThermalEigenstrain
    eigenstrain_name = thermal_eigenstrain
    deformation_gradient_name = thermal_deformation_gradient
    temperature = temperature
    thermal_expansion_coefficients = '12.8e-06 0 0'
    base_name = phase0
  []
  [./strain_phase0]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y'
    base_name = phase0
    eigenstrain_names = thermal_eigenstrain
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
      type = ComputePlaneFiniteStrain
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
    [./disp_x]
      type = StressDivergenceTensors
      variable = disp_x
      eigenstrain_names = thermal_eigenstrain
      component = 0
    [../]
    [./disp_y]
      type = StressDivergenceTensors
      variable = disp_y
      eigenstrain_names = thermal_eigenstrain
        component = 1
    [../]

      [./solid_z]
        type = WeakPlaneStress
        variable = strain_zz
        eigenstrain_names = thermal_eigenstrain
      [../]
    
[]

[Postprocessors]
  [./vonmises]
    type = ElementAverageValue
    variable = vonmises
  [../]
    [./stress_xx]
      type = ElementAverageValue
      variable = stress_xx
    [../]
    [./stress_yy]
        type = ElementAverageValue
        variable = stress_yy
    [../]
    [./stress_zz]
        type = ElementAverageValue
        variable = stress_zz
    [../]
    [./stress_xy]
      type = ElementAverageValue
      variable = stress_xy
    [../]
    [./strain_xx]
      type = ElementAverageValue
      variable = strain_xx
    [../]
    [./strain_xy]
      type = ElementAverageValue
      variable = strain_xy
    [../]
    [./strain_yy]
      type = ElementAverageValue
      variable = strain_yy
    [../]
  [eth_xx]
    type = ElementAverageValue
    variable = eth_xx
  []
  [eth_yy]
    type = ElementAverageValue
    variable = eth_yy
  []
  [eth_zz]
    type = ElementAverageValue
    variable = eth_zz
  []
  [fth_xx]
    type = ElementAverageValue
    variable = fth_xx
  []
  [fth_yy]
    type = ElementAverageValue
    variable = fth_yy
  []
  [fth_zz]
    type = ElementAverageValue
    variable = fth_zz
  []
  [temperature]
    type = ElementAverageValue
    variable = temperature
  []
  [./react_z]
    type = MaterialTensorIntegral
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
  [../]
  [./min_strain_zz]
    type = NodalExtremeValue
    variable = strain_zz
    value_type = min
  [../]
  [./max_strain_zz]
    type = NodalExtremeValue
    variable = strain_zz
    value_type = max
  [../]
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
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
  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-9

  end_time = 50

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 5e-8
    cutback_factor = 0.75
    growth_factor = 1.2
    optimal_iterations = 20
  [../]

[]

[Outputs]
  csv = true
  exodus = true
[]

[Debug]
  show_material_props = true
[]
