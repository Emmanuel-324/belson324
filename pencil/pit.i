#
# KKS model #eta=1 is liquid and eta=0 is solid
#

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    elem_type = QUAD4
    nx = 160
    ny = 80
    xmin = 0
    xmax = 80
    ymin = 0
    ymax = 40
  []
  [toppit]
    input = gmg
    type = SubdomainBoundingBoxGenerator
    bottom_left = '38.5 38 0'
    top_right = '41.5 40 0'
    block_id = 1
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = toppit
    primary_block = '1'
    paired_block = '0'
    new_boundary = 'interface'
  []
  [./break_boundary]
    input = interface
    type = BreakBoundaryOnSubdomainGenerator
  []
[]

[Variables]
  # order parameter
  [eta]
    order = FIRST
    family = LAGRANGE
  []

  # solute concentration
  [c]
    order = FIRST
    family = LAGRANGE
  []

  # solute phase concentration (matrix)
  [cs]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1.000
  []
  # solute phase concentration (precipitate)
  [cl]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0356
  []
[]



[ICs]
  [eta_ic]
    variable = eta
    type = SmoothCircleIC
    x1 = 40
    y1 = 40
    radius = 1.5
    invalue = 1
    outvalue = 0
  []
  [c_ic]
    variable = c
    type = SmoothCircleIC
    x1 = 40
    y1 = 40
    radius = 1.5
    invalue = 0
    outvalue = 1
  []
[]

[BCs]
  [all]
    type =  NeumannBC
    variable = 'c'
    boundary = 'top_to_0 left right bottom'
    value = 0.0
  []	
  [toppit]
    type =  DirichletBC
    variable = 'c'
    boundary = 'top_to_1'
    value = 0.0
  []	
  [bottom]
    type =  DirichletBC
    variable = 'eta'
    boundary = 'bottom'
    value = 0.0
  []	
  [toppiteta]
    type =  DirichletBC
    variable = 'eta'
    boundary = 'top_to_1'
    value = 1.0
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
[]

[Bounds]
  [eta_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta
    bound_type = upper
    bound_value = 1
  []
  [eta_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta
    bound_type = lower
    bound_value = 0
  []
  [c_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = c
    bound_type = lower
    bound_value = 0
  []
  [cl_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = cl
    bound_type = lower
    bound_value = 0
  []
  [cs_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = cs
    bound_type = lower
    bound_value = 0
  []
  [cl_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = cl
    bound_type = upper
    bound_value = 0.1
  []
  [cs_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = cs
    bound_type = upper
    bound_value = 1
  []

[]


[Materials]
  # Chemical free energy of the matrix
  [fs]
    type = DerivativeParsedMaterial
    property_name = fs
    coupled_variables = 'cs'
    expression = '1.61*(cs-1.0)^2'
  []

  # Free energy of the precipitate phase
  [fl]
    type = DerivativeParsedMaterial
    property_name = fl
    coupled_variables = 'cl'
    expression = '1.61*(cl-0.0356)^2'
  []

  # h(eta)
  [h_eta]
    type = SwitchingFunctionMaterial
    h_order = SIMPLE
    eta = eta
  []

  # g(eta)
  [g_eta]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta
  []

  # constant properties
  [constants]
    type = GenericConstantMaterial
    prop_names  =  'L     kappa   Dl   Ds'
    prop_values = '1.957  0.0833    1   1e-4'  
  []

  # Coefficients for diffusion equation
  [D]
    type = DerivativeParsedMaterial
    material_property_names = 'Dl Ds h'
    coupled_variables = 'eta'
    expression = Dl*h+Ds*(1-h)
    property_name = D
  []

  # Coefficients for diffusion equation
  [Dhs]
    type = DerivativeParsedMaterial
    material_property_names = 'D h'
    expression = D*(1-h)
    property_name = Dhs
  []
  [Dhl]
    type = DerivativeParsedMaterial
    material_property_names = 'D h'
    expression = D*h
    property_name = Dhl
  []

[]

[Kernels]

  # enforce c = (1-h(eta))*cs + h(eta)*cl
  [PhaseConc]
    type = KKSPhaseConcentration
    variable = cl
    ca       = cs
    c        = c
    eta      = eta
  []

  # enforce pointwise equality of chemical potentials
  [ChemPotVacancies]
    type = KKSPhaseChemicalPotential
    variable = cs
    cb       = cl
    fa_name  = fs
    fb_name  = fl
  []

  #Kernels for diffusion equation
  [diff_time]
    type = TimeDerivative
    variable = c
  []
  [diff_cl]
    type = MatDiffusion
    variable = c
    diffusivity = Dl
    v = cl
  []
  [diff_cs]
    type = MatDiffusion
    variable = c
    diffusivity = Ds
    v = cs
  []

  #
  # Allen-Cahn Equation
  #
  [ACBulkF]
    type = KKSACBulkF
    variable = eta
    fa_name  = fs
    fb_name  = fl
    w        = 1.0
    coupled_variables = 'cl cs'
  [../]
  [./ACBulkC]
    type = KKSACBulkC
    variable = eta
    ca       = cs
    cb       = cl
    fa_name  = fs
  []
  [ACInterface]
    type = ACInterface
    variable = eta
    kappa_name = kappa
  []

  [detadt]
    type = TimeDerivative
    variable = eta
  []
[]

[VectorPostprocessors]
  [eta1]
    type = LineValueSampler
    start_point = '40 40 0'
    end_point = '40 0 0'
    num_points = 400
    sort_by = y
    variable = eta
  []
  [c]
    type = LineValueSampler
    start_point = '40 40 0'
    end_point = '40 0 0'
    num_points = 400
    sort_by = y
    variable = c
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  petsc_options_value = 'asm       lu            nonzero'

  l_max_its = 50
  nl_max_its = 25
  l_tol = 1.0e-3
  nl_rel_tol = 1.0e-7
  nl_abs_tol = 1.0e-9

  end_time = 60000

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.0005
    cutback_factor = 0.75
    growth_factor = 1.2
    optimal_iterations = 20
  []

  [Adaptivity]
    initial_adaptivity = 0
    refine_fraction = 0.7
    coarsen_fraction = 0.1
    max_h_level = 2
  []
[]
#
# Precondition using handcoded off-diagonal terms
#
[Preconditioning]
  [full]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  [dofs]
    type = NumDOFs
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
