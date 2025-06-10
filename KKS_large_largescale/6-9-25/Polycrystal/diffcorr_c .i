[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [./gen]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    nx = 80
    ny = 40
    xmin = 0
    xmax = 40
    ymin = 0
    ymax = 20
  [../]
  [./grain1]
    input = gen
    type = SubdomainBoundingBoxGenerator
    block_id = 0
    bottom_left = '0 0 0'
    top_right = '20 20 0'
  [../]
  [./grain2]
    input = grain1
    type = SubdomainBoundingBoxGenerator
    block_id = 1
    bottom_left = '20 0 0'
    top_right = '40 20 0'
  [../]
  [./toppit]
    input = grain2
    type = SubdomainBoundingBoxGenerator
    bottom_left = '18.5 19.75 0'
    top_right = '21.5 20 0'
    block_id = 2
  [../]
  [./interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = toppit
    primary_block = '1'
    paired_block = '0'
    new_boundary = 'interface'
  [../]
  [./break_boundary]
    input = interface
    type = BreakBoundaryOnSubdomainGenerator
  [../]
[]

[Variables]
  # order parameter
  [./eta]
    order = FIRST
    family = LAGRANGE
  [../]

  # solute concentration
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]

  # solute phase concentration (matrix)
  [./cs]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1.000
  [../]
  # solute phase concentration (precipitate)
  [./cl]
    order = FIRST
    family = LAGRANGE
#    initial_condition = 0.0356
  [../]
[]

[ICs]
  [./eta_ic]
    variable = eta
    type = SmoothCircleIC
    x1 = 20
    y1 = 20
    radius = 1.5
#    int_width = 0.05
    invalue = 1
    outvalue = 0
  [../]
  [./c_ic]
    variable = c
    type = SmoothCircleIC
    x1 = 20
    y1 = 20
    radius = 1.5
#    int_width = 0.05
    invalue = 0
    outvalue = 1
  [../]
[]



[BCs]
  [./all]
    type =  NeumannBC
    variable = 'c'
    boundary = 'top_to_0 top_to_1 left right bottom'
    value = 0.0
  [../]	
  [./toppit]
    type =  DirichletBC
    variable = 'c'
    boundary = 'top_to_2'
    value = 0.0
  [../]	
  [./bottom]
    type =  DirichletBC
    variable = 'eta'
    boundary = 'bottom'
    value = 0.0
  [../]	
  [./toppiteta]
    type =  DirichletBC
    variable = 'eta'
    boundary = 'top_to_2'
    value = 1.0
  [../]	
  [./stressfree_boundary]
     #Applies the pressure
     type = Pressure
     boundary = top
     factor = 0.0 
     variable = disp_y
  [../]
  [./dispy]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./dispx]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./tdisp]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = tdisp
  [../]
[]

[UserObjects]
  [./prop_read]
    type = PropertyReadFile
    prop_file_name = 'euler_ang_file.txt'
    # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 3
    read_type = block
    nblock= 3
  [../]
[]


[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]
  [./accslip]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss1]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]


[Bounds]
  [./eta_upper_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = eta
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_lower_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = eta
    bound_type = lower
    bound_value = 0
  [../]
  [./c_lower_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = c
    bound_type = lower
    bound_value = 0
  [../]
  [./cl_lower_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = cl
    bound_type = lower
    bound_value = 0
  [../]
  [./cs_lower_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = cs
    bound_type = lower
    bound_value = 0
  [../]
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    value = 2.9400e-07*t  #2.5e-8(strain_rate)*20(x_length)*tc(0.0294)
  [../]
[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
  add_variables = true
  strain_base_name = uncracked
[]

[AuxKernels]
  [./accslip]
    type = MaterialRealAux
    variable = accslip
    property = acc_slip
    execute_on = timestep_end
#   block = 0
  [../]
  [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    scalar_type = VonMisesStress
    execute_on = timestep_end
#   block = 0
  [../]
  [./fp_xx]
    type = RankTwoAux
    variable = fp_xx
    rank_two_tensor = fp
    index_j = 0
    index_i = 0
    execute_on = timestep_end
#    block = 0
  [../]
  [./e_xx]
    type = RankTwoAux
    variable = e_xx
    rank_two_tensor = lage
    index_j = 0
    index_i = 0
    execute_on = timestep_end
#    block = 0
  [../]
  [./gss1]
    type = MaterialStdVectorAux
    variable = gss1
    property = gss
    index = 0
    execute_on = timestep_end
#    block = 0
  [../]
[]


[Materials]
  # Chemical free energy of the matrix
  [./fs]
    type = DerivativeParsedMaterial
    f_name = fs
    args = 'cs'
    function = '1.61*(cs-1.0)^2'
  [../]

  # Free energy of the precipitate phase
  [./fl]
    type = DerivativeParsedMaterial
    f_name = fl
    args = 'cl'
    function = '1.61*(cl-0.0356)^2'
  [../]

  # h(eta)
  [./h_eta]
    type = SwitchingFunctionMaterial
    h_order = SIMPLE
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
    prop_names  =  'L0        kappa   Dl   Ds'
    prop_values = '1.957e-4  0.0833    1   1e-4'  
  [../]

  [./reaction_mob]
    type = ParsedMaterial
    f_name = L
    args = 'accslip'
    material_property_names = 'L0'
    function = 'L0*(1+accslip/0.001)'
  [../]

  # Coefficients for diffusion equation
  [./D]
    type = DerivativeParsedMaterial
    material_property_names = 'Dl Ds h'
    args = 'eta'
    function = Dl*h+Ds*(1-h)
    f_name = D
  [../]

  # Coefficients for diffusion equation
  [./Dhs]
    type = DerivativeParsedMaterial
    material_property_names = 'D h'
    function = D*(1-h)
    f_name = Dhs
  [../]
  [./Dhl]
    type = DerivativeParsedMaterial
    material_property_names = 'D h'
    function = D*h
    f_name = Dhl
  [../]

  [./crysp]
    type = FiniteStrainCPSlipRateRes
    slip_sys_file_name = input_slip_sys.txt
    nss = 12
    num_slip_sys_flowrate_props = 2 #Number of properties in a slip system
    flowprops = '1 4 1e-6 0.05 5 8 1e-6 0.05 9 12 1e-6 0.05'
    hprops = '1.0 215 78.6 200.0 1.75'
    gprops = '1 4 78.6 5 8 78.6 9 12 78.6'
    slip_incr_tol = 1
    maximum_substep_iteration = 12
    rtol = 1e-9
    abs_tol = 1e-10
    base_name = uncracked
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP
    C_ijkl = '2.04e5 1.37e5 1.37e5 2.044e5 1.37e5 2.04e5 1.26e5 1.26e5 1.26e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
    output_properties = 'Euler_angles'
    outputs = exodus
    base_name = uncracked
  [../]
  [./cracked_stress]
    type = ComputeCrackedStress
    c = eta
    F_name = E_el
    use_current_history_variable = false
    uncracked_base_name = uncracked
    finite_strain_model = false
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'eta'
#    function = '1'
    function = '(1-eta^3*(6*eta^2-15*eta+10))*(1.0 - c) + c'
    constant_names       = 'c'
    constant_expressions = '1e-6'
    derivative_order = 2
  [../]

[]


[Kernels]
  # enforce c = (1-h(eta))*cs + h(eta)*cl
  [./PhaseConc]
    type = KKSPhaseConcentration
    variable = cl
    ca       = cs
    c        = c
    eta      = eta
  [../]

  # enforce pointwise equality of chemical potentials
  [./ChemPotVacancies]
    type = KKSPhaseChemicalPotential
    variable = cs
    cb       = cl
    fa_name  = fs
    fb_name  = fl
  [../]

  #Kernels for diffusion equation
  [./diff_time]
    type = TimeDerivative
    variable = c
  [../]
  [./diff_cl]
    type = MatDiffusion
    variable = c
    diffusivity = Dl
    v = cl
  [../]
  [./diff_cs]
    type = MatDiffusion
    variable = c
    diffusivity = Ds
    v = cs
  [../]

  #
  # Allen-Cahn Equation
  #
  [./ACBulkF]
    type = KKSACBulkF
    variable = eta
    fa_name  = fs
    fb_name  = fl
    w        = 1.0
    args = 'cl cs'
  [../]
  [./ACBulkC]
    type = KKSACBulkC
    variable = eta
    ca       = cs
    cb       = cl
    fa_name  = fs
  [../]
  [./ACInterface]
    type = ACInterface
    variable = eta
    kappa_name = kappa
  [../]

  [./detadt]
    type = TimeDerivative
    variable = eta
  [../]
[]

[VectorPostprocessors]
  [./eta1]
    type = LineValueSampler
    start_point = '20 20 0'
    end_point = '20 0 0'
    num_points = 200
    sort_by = y
    variable = eta
  [../]
  [./c]
    type = LineValueSampler
    start_point = '20 20 0'
    end_point = '20 0 0'
    num_points = 200
    sort_by = y
    variable = c
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu            superlu_dist          vinewtonrsls'  

  l_max_its = 20
  nl_max_its = 10
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-9

  dtmax = 350
  end_time = 100000

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.05
    cutback_factor = 0.75
    growth_factor = 1.2
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
    refine_fraction = 0.6 # Fraction of high error that will be refined
    coarsen_fraction = 0.1 # Fraction of low error that will coarsened
    max_h_level = 2 # Max number of refinements used, starting from initial mesh (befo# uniform refinement)
  [../]
[]

# Precondition using handcoded off-diagonal terms

[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./dofs]
    type = NumDOFs
  [../]
[]

[Outputs]
  exodus = true
  interval = 2
  [./console]
     type = Console
     interval = 2
  [../]
#  [./csv]
#     type = CSV
#     interval = 2
#  [../]
[]

