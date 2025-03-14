[Mesh]
    [line_1D]
      type = GeneratedMeshGenerator
      dim = 1
      nx = 10
      xmin = 0
      xmax = 1
    []
  []
  
  [GlobalParams]
    displacements = 'disp_x'
  []
  
  [Variables]
    [disp_x]
      family = LAGRANGE
      order = FIRST
    []
  []
  
  [Functions]
    [compress_func]
      type = ParsedFunction
      expression = '50e6 * t'
    []
  []
  
  [Physics]
    [SolidMechanics]
      [QuasiStatic]
        [all]
          strain = SMALL
          incremental = true
          add_variables = true
          generate_output = 'stress_xx strain_xx vonmises_stress PLASTIC_STRAIN_XX effective_plastic_strain'
          material_output_family = MONOMIAL
          material_output_order = FIRST
        []
      []
    []
  []
  
  [BCs]
    [left_fix]
      type = DirichletBC
      variable = disp_x
      boundary = 'left'
      value = 0.0
    []
    [right_pressure]
      type = Pressure
      variable = disp_x
      boundary = right
      function = 'compress_func'
    []
  []
  
  [ICs]
    [ic_disp_x]
      type = ConstantIC
      variable = disp_x
      value = 0.0
    []
  []
  
  [UserObjects]
    [cohesion_obj]
      type = SolidMechanicsHardeningConstant
      value = 10e6
    []
    [phi_obj]
      type = SolidMechanicsHardeningConstant
      value = '35'
      convert_to_radians = true
    []
    [psi_obj]
      type = SolidMechanicsHardeningConstant
      value = '5'
      convert_to_radians = true
    []
    [dp_model]
      type = SolidMechanicsPlasticDruckerPrager
      mc_cohesion = 'cohesion_obj'
      mc_friction_angle = 'phi_obj'
      mc_dilation_angle = 'psi_obj'
      yield_function_tolerance = '1e-10'
      internal_constraint_tolerance = '1e-10'
    []
  []
  
  [Materials]
    [elasticity]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 15e9
      poissons_ratio = 0.2
      block = 0
    []
    [drucker_prager]
      type = ComputeMultiPlasticityStress
      ep_plastic_tolerance = 1e-5
      plastic_models = dp_model # Defined in UserObjects, Mohr-Coulomb yield surface
    []
  []
  
  [Executioner]
    type = Transient
    solve_type = 'PJFNK'
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
    petsc_options_value = 'lu            superlu_dist          vinewtonrsls' 
    end_time = 10
    dt = 5e-4
    # nl_abs_tol = 1e-30
    # nl_rel_tol = 1e-50
    # nl_max_its = 200
  []
  [Postprocessors]
    
     [./stress_xx]
       type = ElementAverageValue
       variable = stress_xx
     [../]
    [./strain_xx]
        type = ElementAverageValue
        variable = strain_xx
    [../]
[]
  [Outputs]
    [Exodus]
      type = Exodus
    []
    [./table]
        type = CSV
        execute_on = timestep_end
    [../]
  []