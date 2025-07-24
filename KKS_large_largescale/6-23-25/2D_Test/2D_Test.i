[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [copper]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    nx = 175
    ny = 175
    xmin = 0
    xmax = 350
    ymin = 0
    ymax = 350
  []
  [copper_id]
    type = SubdomainIDGenerator
    input = copper
    subdomain_id = 0
  []

  [brass]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    nx = 175
    ny = 175
    xmin = 350     # Start exactly where copper ends
    xmax = 700
    ymin = 0
    ymax = 350
  []
  [brass_id]
    type = SubdomainIDGenerator
    input = brass
    subdomain_id = 1
  []

  [sticher]
    type = StitchedMeshGenerator
    inputs = 'copper_id brass_id'
    stitch_boundaries_pairs = 'right left'
    prevent_boundary_ids_overlap = false
  []
[]


[AuxVariables]
  [pk2]
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
  [copper_gss]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  []
  [copper_slip_increment]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  []
  [brass_gss]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
  [brass_slip_increment]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [copper]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = stress_xx
    block = 0
    base_name = copper
  []
  [brass]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = stress_xx
    block = 1
    base_name = brass
  []
[]

[AuxKernels]
  [pk2]
    type = RankTwoAux
    variable = pk2
    rank_two_tensor = second_piola_kirchhoff_stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []
  [fp_xx]
    type = RankTwoAux
    variable = fp_xx
    rank_two_tensor = plastic_deformation_gradient
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []
  [e_xx]
    type = RankTwoAux
    variable = e_xx
    rank_two_tensor = total_lagrangian_strain
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []
  [gss_copper]
    type = MaterialStdVectorAux
    variable = copper_gss
    property = copper_slip_resistance
    index = 0
    block = 0
    execute_on = timestep_end
  []
  [slip_inc_copper]
    type = MaterialStdVectorAux
    variable = copper_slip_increment
    property = copper_slip_increment
    index = 0
    block = 0
    execute_on = timestep_end
  []
  [gss_brass]
    type = MaterialStdVectorAux
    variable = brass_gss
    property = brass_slip_resistance
    index = 0
    block = 1
    execute_on = timestep_end
  []
  [slip_inc_brass]
    type = MaterialStdVectorAux
    variable = brass_slip_increment
    property = brass_slip_increment
    index = 0
    block = 1
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
    function = '0.01*t'
  []
[]

[Materials]
  [elasticity_tensor_copper]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
    base_name = copper
    block = 0
  []
  [stress_copper]
    type = ComputeMultipleCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl_copper'
    tan_mod_type = exact
    base_name = copper
    block = 0
  []
  [trial_xtalpl_copper]
    type = CrystalPlasticityKalidindiUpdate
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
    base_name = copper
    block = 0
  []
  [elasticity_tensor_brass]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
    euler_angle_1 = 0.0
    euler_angle_2 = 45.0
    euler_angle_3 = 0.9
    base_name = brass
    block = 1
  []
  [stress_brass]
    type = ComputeMultipleCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl_brass'
    tan_mod_type = exact
    base_name = brass
    block = 1
  []
  [trial_xtalpl_brass]
    type = CrystalPlasticityKalidindiUpdate
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
    base_name = brass
    block = 1
  []
[]

[Postprocessors]
  [copper_stress_xx]
    type = ElementAverageValue
    variable = copper_stress_xx
    block = 0
  []
  [brass_stress_xx]
    type = ElementAverageValue
    variable = brass_stress_xx
    block = 1
  []
  [pk2]
    type = ElementAverageValue
    variable = pk2
  []
  [fp_xx]
    type = ElementAverageValue
    variable = fp_xx
  []
  [e_xx]
    type = ElementAverageValue
    variable = e_xx
  []
  [copper_gss]
    type = ElementAverageValue
    variable = copper_gss
    block = 0
  []
  [copper_slip_increment]
    type = ElementAverageValue
    variable = copper_slip_increment
    block = 0
  []
  [brass_gss]
    type = ElementAverageValue
    variable = brass_gss
    block = 1
  []
  [brass_slip_increment]
    type = ElementAverageValue
    variable = brass_slip_increment
    block = 1
  []
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

  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = ' asm      2              lu            gmres     200'
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-10
  nl_abs_step_tol = 1e-10

  dt = 0.05
  dtmin = 0.01
  dtmax = 10.0
  num_steps = 10
[]

[Outputs]
  csv = true
  exodus = true
  perf_graph = true
[]
