[GlobalParams]
  dim = 2
[../]
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
 

 
 [AuxVariables]
   [./vonmises]
     order = CONSTANT
     family = MONOMIAL
   [../]
   [./stress_xx]
     order = CONSTANT
     family = MONOMIAL
   [../]
    [./stress_zz]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [temperature]
      order = FIRST
      family = LAGRANGE
    []
    
    [feig_xx]
      order = CONSTANT
      family = MONOMIAL
    []
    [feig_yy]
      order = CONSTANT
      family = MONOMIAL
    []
    [f1_xx]
      order = CONSTANT
      family = MONOMIAL
    []
    [f1_yy]
      order = CONSTANT
      family = MONOMIAL
    []
 []

 [AuxKernels]
   [vonmises]
     type = RankTwoScalarAux
     rank_two_tensor = stress
     variable = vonmises
     scalar_type = VonMisesStress
     execute_on = timestep_end
   [../]
   [stress_xx]
     type = RankTwoAux
     rank_two_tensor = stress
     variable = stress_xx
     index_j = 0
     index_i = 0
     execute_on = timestep_end 
   [../]
    [stress_zz]
      type = RankTwoAux
      rank_two_tensor = stress
      variable = stress_zz
      index_j = 2
      index_i = 2
      execute_on = timestep_end 
    [../]
     
    
    [feig_xx]
      type = RankTwoAux
      variable = feig_xx
      rank_two_tensor = eigenstrain_deformation_gradient
      index_j = 0
      index_i = 0
      execute_on = timestep_end
    []
    [feig_yy]
      type = RankTwoAux
      variable = feig_yy
      rank_two_tensor = eigenstrain_deformation_gradient
      index_j = 1
      index_i = 1
      execute_on = timestep_end
    []
    [temperature]
      type = FunctionAux
      variable = temperature
      function = '300'
      execute_on = 'timestep_begin'
    []
    [f1_xx]
      type = RankTwoAux
      variable = f1_xx
      rank_two_tensor = phase0_deformation_gradient_1
      index_j = 0
      index_i = 0
      execute_on = timestep_end
    []
    [f1_yy]
      type = RankTwoAux
      variable = f1_yy
      rank_two_tensor = phase0_deformation_gradient_1
      index_j = 1
      index_i = 1
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
  
  [stress_free_bc]
    type = Pressure
    variable = disp_x
    boundary = right
    displacements = 'disp_x disp_y'
    factor = 0
  []

  [stress_free_bc2]
    type = Pressure
    variable = disp_y
    boundary = top
    displacements = 'disp_x disp_y'
    factor = 0
  []
 
  
 []
 
 [Materials]
 
  [elasticity_tensor_phase0]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.906e5
    poissons_ratio = 0.3
    base_name = phase0
    block = 0
  []
  [stress_phase0]
    type = ComputeMultipleCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl_phase0'
    eigenstrain_names = 'eigenstrain_0'
    tan_mod_type = exact
    rtol = 1e-08
    base_name = phase0
  []
  
  [trial_xtalpl_phase0]
    type = CrystalPlasticityKalidindiUpdate
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
    crystal_lattice_type = FCC
    resistance_tol = 0.01
    r = 1.4             
    h = 6000            
    t_sat = 598.5        
    gss_a = 1.5         
    ao = 0.001           
    xm = 0.017             
    gss_initial = 600
    base_name = phase0
  []
  [eigenstrain_0]
    type = ComputeCrystalPlasticityThermalEigenstrain
    eigenstrain_name = eigenstrain_0
    deformation_gradient_name = deformation_gradient_1
    temperature = temperature
    thermal_expansion_coefficients = '12.8e-06 12.8e-06 12.8e-06'
    base_name = phase0
  []
  
  [./strain_phase0]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y'
    base_name = phase0
  [../]

  
  [elasticity_tensor_phase1]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.25
    base_name = phase1
  []
  [stress_phase1]
    type = ComputeMultipleCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl_phase1'
    tan_mod_type = exact
    rtol = 1e-08
    base_name = phase1
  []
  [trial_xtalpl_phase1]
    type = CrystalPlasticityKalidindiUpdate
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
   [./TensorMechanics]
      displacements = 'disp_x disp_y'
      add_variables = true
      generate_output = 'stress_xx stress_zz'
      strain = FINITE
   []
    
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
    [./stress_zz]
      type = ElementAverageValue
      variable = stress_zz
    [../]
    [f1_xx]
      type = ElementAverageValue
      variable = f1_xx
    []
    [f1_yy]
      type = ElementAverageValue
      variable = f1_yy
    []
    [./feig_xx]
      type = ElementAverageValue
      variable = feig_xx
    [../]
    [./feig_yy]
      type = ElementAverageValue
      variable = feig_yy
    [../]
    [temperature]
      type = ElementAverageValue
      variable = temperature
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
   l_max_its = 50
   nl_max_its = 50
   nl_rel_tol = 1.0e-8
   nl_abs_tol = 1.0e-9
 
   end_time = 100
 
   [./TimeStepper]
     type = IterationAdaptiveDT
     dt = 5e-8
     cutback_factor = 0.75
     growth_factor = 1.2
     optimal_iterations = 50
   [../]
 
 []
 
 [Outputs]
   exodus = true
   csv = true
 []