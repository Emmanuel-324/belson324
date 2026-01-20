[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [file]
     type = FileMeshGenerator
     file = CR_1_out.e
     use_for_exodus_restart = true
   []
 []
 
 [Variables]
   [./disp_x]
   [../]
   [./disp_y]
   [../]
 
   # order parameter 0
   [./eta_m]
     initial_from_file_var = eta_m
   [../]
   # order parameter 1
   [./eta_gp]
     initial_from_file_var = eta_gp
   [../]
    # order parameter 2
    [./eta_gpp1]
      initial_from_file_var = eta_gpp1
    [../]
    [./eta_gpp2]
      initial_from_file_var = eta_gpp2
    [../]  
 []
 
 [AuxVariables]
   [./vonmises]
     order = CONSTANT
     family = MONOMIAL
   [../]
 []
 
 [Physics/SolidMechanics/QuasiStatic/all]
  strain = FINITE
  incremental = true
  add_variables = true
  generate_output = 'stress_xx stress_xy stress_yy stress_zz strain_xx strain_xy strain_yy strain_zz'
 []

 
 [AuxKernels]
   [./vonmises]
     type = RankTwoScalarAux
     rank_two_tensor = stress
     variable = vonmises
     scalar_type = VonMisesStress
     execute_on = timestep_end
   [../]
   
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
    [elasticity_tensor_phasem]
      type = ComputeElasticityTensorCP
      C_ijkl = '1698.5 1054.8 1054.8 1698.5 1054.8 1698.5 818.6 818.6 818.6' #Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of precipitates
      fill_method = symmetric9
      base_name = phasem
    []
    [stress_phasem]
      type = ComputeMultipleCrystalPlasticityStress_abs
      crystal_plasticity_models = 'trial_xtalpl_phasem'
      tan_mod_type = exact
      rtol = 1e-08
      base_name = phasem
    []
    [trial_xtalpl_phasem]
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
      gss_initial = 665.5 
      base_name = phasem
    []
   
    [./strain_phasem]
      type = ComputeFiniteStrain
      displacements = 'disp_x disp_y'
      base_name = phasem
 #    eigenstrain_names = eigenstrain2
    [../]
    
    [elasticity_tensor_phasegp]
      type = ComputeElasticityTensorCP
      C_ijkl = '1516.7 966.2 966.2 1516.7 966.2 1516.7 825.7 825.7 825.7' #Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of precipitates
      fill_method = symmetric9
      base_name = phasegp
    []
    [stress_phasegp]
      type = ComputeMultipleCrystalPlasticityStress_abs
      crystal_plasticity_models = 'trial_xtalpl_phasegp'
      tan_mod_type = exact
      rtol = 1e-08
      base_name = phasegp
    []
    [trial_xtalpl_phasegp]
      type = CrystalPlasticityKalidindiUpdate_slip
      slip_sys_file_name = input_slip_sys.txt
      crystal_lattice_type = FCC
      resistance_tol = 0.01
      number_slip_systems = 12
      slip_system_modes = 3
      number_slip_systems_per_mode = '4 4 4'
      lattice_friction_per_mode = '1.872  3.483 4.065'
      r = 1.0             
      h = 37.45  # Aitor Cruzado et. al  DOI: 10.1007/978-3-030-40562-5_5     
      t_sat = 3.735  # Aitor Cruzado et. al  DOI: 10.1007/978-3-030-40562-5_5      
      gss_a = 1.5         
      ao = 0.001           
      xm = 0.017 # Aitor Cruzado et. al  DOI: 10.1007/978-3-030-40562-5_5            
      base_name = phasegp
    []
    [./strain_phasegp]
      type = ComputeFiniteStrain
      displacements = 'disp_x disp_y'
      base_name = phasegp
    [../]


    [elasticity_tensor_phasegpp1]
      type = ComputeElasticityTensorCP
      C_ijkl = '1814.2 1167.2  1167.2  1814.2 1167.2 1814.2 712.5 712.5 712.5' #Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of precipitates
      fill_method = symmetric9
      base_name = phasegpp1
    []
    [stress_phasegpp1]
      type = ComputeMultipleCrystalPlasticityStress_abs
      crystal_plasticity_models = 'trial_xtalpl_phasegpp1'
      tan_mod_type = exact
      rtol = 1e-08
      base_name = phasegpp1
    []
    [trial_xtalpl_phasegpp1]
      type = CrystalPlasticityKalidindiUpdate_slip
      slip_sys_file_name = input_slip_sys.txt
      crystal_lattice_type = FCC
      resistance_tol = 0.01
      number_slip_systems = 12
      slip_system_modes = 3
      number_slip_systems_per_mode = '4 4 4'
      lattice_friction_per_mode = '1.872  3.483 4.065'
      r = 1.0             
      h = 37.45  # Aitor Cruzado et. al  DOI: 10.1007/978-3-030-40562-5_5     
      t_sat = 3.735  # Aitor Cruzado et. al  DOI: 10.1007/978-3-030-40562-5_5      
      gss_a = 1.5         
      ao = 0.001           
      xm = 0.017 # Aitor Cruzado et. al  DOI: 10.1007/978-3-030-40562-5_5            
      base_name = phasegpp1
    []
    [./strain_phasegpp1]
      type = ComputeFiniteStrain
      displacements = 'disp_x disp_y'
      base_name = phasegpp1
    [../]

    [elasticity_tensor_phasegpp2]
      type = ComputeElasticityTensorCP
      C_ijkl = '1814.2 1167.2  1167.2  1814.2 1167.2 1814.2 712.5 712.5 712.5' #Ghorbanpour, S., et al., A crystal plasticity model incorporating the effects of precipitates
      fill_method = symmetric9
      base_name = phasegpp2
    []
    [stress_phasegpp2]
      type = ComputeMultipleCrystalPlasticityStress_abs
      crystal_plasticity_models = 'trial_xtalpl_phasegpp2'
      tan_mod_type = exact
      rtol = 1e-08
      base_name = phasegpp2
    []
    [trial_xtalpl_phasegpp2]
     type = CrystalPlasticityKalidindiUpdate_slip
      slip_sys_file_name = input_slip_sys.txt
      crystal_lattice_type = FCC
      resistance_tol = 0.01
      number_slip_systems = 12
      slip_system_modes = 3
      number_slip_systems_per_mode = '4 4 4'
      lattice_friction_per_mode = '1.872  3.483 4.065'
      r = 1.0             
      h = 37.45  # Aitor Cruzado et. al  DOI: 10.1007/978-3-030-40562-5_5     
      t_sat = 3.735  # Aitor Cruzado et. al  DOI: 10.1007/978-3-030-40562-5_5      
      gss_a = 1.5         
      ao = 0.001           
      xm = 0.017 # Aitor Cruzado et. al  DOI: 10.1007/978-3-030-40562-5_5            
      base_name = phasegpp2
    []
    [./strain_phasegpp2]
      type = ComputeFiniteStrain
      displacements = 'disp_x disp_y'
      base_name = phasegpp2
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
  
   # Generate the global stress from the phase stresses
    [./global_stress]
      type = MultiPhaseStressMaterial
      phase_base = 'phasegp phasegpp1 phasegpp2 phasem'
      h          = 'hgp     hgpp1   hgpp2   hm'
    [../]
 
[]
 
 [Kernels]
   [./etam_dt]
     type = TimeDerivative
     variable = eta_m
   [../]
   [./etagp_dt]
     type = TimeDerivative
     variable = eta_gp
   [../]
    [./etagpp1_dt]
     type = TimeDerivative
     variable = eta_gpp1
   [../]
    [./etagpp2_dt]
     type = TimeDerivative
     variable = eta_gpp2
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
    [./strain_zz]
      type = ElementAverageValue
      variable = strain_zz
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
 
 [Outputs]
   exodus = true
   csv = true
 []