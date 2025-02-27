#
# KKS model #eta=1 is liquid and eta=0 is solid
#

[Mesh]
  type = GeneratedMesh
  dim = 2
  elem_type = QUAD4
  nx = 200
  ny = 10
  nz = 0
  xmin = 0
  xmax = 100
  ymin = 0
  ymax = 5
  zmin = 0
  zmax = 0
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
    initial_condition = 0.0356
  [../]
[]



[ICs]
  [./eta_ic]
    variable = eta
    type = BoundingBoxIC
    x1 = 0
    x2 = 5
    y1 = 0
    y2 = 5
    inside = 1
    outside = 0
  [../]
  [./c_ic]
    variable = c
    type = BoundingBoxIC
    x1 = 0
    x2 = 5
    y1 = 0
    y2 = 5
    inside = 0.0
    outside = 1.0
  [../]
[]

[BCs]
  [./all]
    type =  NeumannBC
    variable = 'c'
    boundary = 'right top bottom'
    value = 0.0
  [../]	
  [./left]
    type =  DirichletBC
    variable = 'c'
    boundary = 'left'
    value = 0.0
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
  [./cl_upper_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = cl
    bound_type = upper
    bound_value = 0.1
  [../]
  [./cs_upper_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = cs
    bound_type = upper
    bound_value = 1
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
    prop_names  =  'L     kappa   Dl   Ds'
    prop_values = '1.957  0.0833    1   1e-4'  
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
  [eta1]
    type = LineValueSampler
    start_point = '0 2.5 0'
    end_point = '100 2.5 0'
    num_points = 1000
    sort_by = x
    variable = eta
  []
  [c]
    type = LineValueSampler
    start_point = '0 2.5 0'
    end_point = '100 2.5 0'
    num_points = 1000
    sort_by = x
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

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.0005
    cutback_factor = 0.75
    growth_factor = 1.2
    optimal_iterations = 20
  [../]

  [./Adaptivity]
    initial_adaptivity = 0
    refine_fraction = 0.7
    coarsen_fraction = 0.1
    max_h_level = 2
  [../]
[]


#
# Precondition using handcoded off-diagonal terms
#
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
  [./csv]
     type = CSV
     interval = 2
  [../]
[]
