#
# KKS ternary (3 chemical component) system example in the split form
# We track c1 and c2 only, since c1 + c2 + c3 = 1
#

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 150
  ny = 15
  nz = 0
  xmin = -25
  xmax = 25
  ymin = -2.5
  ymax = 2.5
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]


[Variables]
  # order parameter
  [./eta1]
    order = FIRST
    family = LAGRANGE
  [../]
  [./eta2]
    order = FIRST
    family = LAGRANGE
  [../]

  # solute 1 concentration
  [./c1]
    order = FIRST
    family = LAGRANGE
  [../]

  # solute 2 concentration
  [./c2]
    order = FIRST
    family = LAGRANGE
  [../]


  # chemical potential solute 1
  [./w1]
    order = FIRST
    family = LAGRANGE
  [../]

  # chemical potential solute 2
  [./w2]
    order = FIRST
    family = LAGRANGE
  [../]


  # Liquid phase solute 1 concentration
  [./c1l]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.1
  [../]

  # Liquid phase solute 2 concentration
  [./c2l]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.05
  [../]

  # Solid phase solute 1 concentration
  [./c1s]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.8
  [../]

  # Solid phase solute 2 concentration
  [./c2s]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.1
  [../]

[]

[Functions]
  [./ic_func_eta1]
    type = ParsedFunction
    expression = '0.5*(1.0-tanh((x)/sqrt(2.0)))'
  [../]
  [./ic_func_c1]
    type = ParsedFunction
    expression = '0.8*(0.5*(1.0-tanh(x/sqrt(2.0))))^3*(6*(0.5*(1.0-tanh(x/sqrt(2.0))))^2-15*(0.5*(1.0-tanh(x/sqrt(2.0))))+10)+0.1*(1-(0.5*(1.0-tanh(x/sqrt(2.0))))^3*(6*(0.5*(1.0-tanh(x/sqrt(2.0))))^2-15*(0.5*(1.0-tanh(x/sqrt(2.0))))+10))'
  [../]
  [./ic_func_c2]
    type = ParsedFunction
    expression = '0.1*(0.5*(1.0-tanh(x/sqrt(2.0))))^3*(6*(0.5*(1.0-tanh(x/sqrt(2.0))))^2-15*(0.5*(1.0-tanh(x/sqrt(2.0))))+10)+0.05*(1-(0.5*(1.0-tanh(x/sqrt(2.0))))^3*(6*(0.5*(1.0-tanh(x/sqrt(2.0))))^2-15*(0.5*(1.0-tanh(x/sqrt(2.0))))+10))'
  [../]

[]


[ICs]
  [./eta1]
    variable = eta1
    type = FunctionIC
    function = ic_func_eta1
  [../]
  [./eta2]
    variable = eta2
    type = FunctionIC
    function = ic_func_eta1
  [../]
  [./c1]
    variable = c1
    type = FunctionIC
    function = ic_func_c1
  [../]
  [./c2]
    variable = c2
    type = FunctionIC
    function = ic_func_c2
  [../]

[]

[Materials]
  # Free energy of the liquid
  [./fl]
    type = DerivativeParsedMaterial
    property_name = fl
    coupled_variables = 'c1l c2l'
    expression = '(0.1-c1l)^2+(0.05-c2l)^2'
  [../]

  # Free energy of the solid
  [./fs]
    type = DerivativeParsedMaterial
    property_name = fs
    coupled_variables = 'c1s c2s'
    expression = '(0.8-c1s)^2+(0.1-c2s)^2'
  [../]

  # h(eta1)
  [./h_eta1]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = h1
    all_etas = 'eta1 eta2'
    phase_etas = eta1
  [../]

  # g(eta1)
  [./g_eta1]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta1
    function_name = g1
  [../]

 # h(eta2)
 [./h_eta2]
  type = SwitchingFunctionMultiPhaseMaterial
  h_name = h2
  all_etas = 'eta1 eta2'
  phase_etas = eta2
[../]

# g(eta2)
[./g_eta2]
  type = BarrierFunctionMaterial
  g_order = SIMPLE
  eta = eta2
  function_name = g2
[../]

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'M   L   eps_sq'
    prop_values = '0.7 0.7 1.0  '
  [../]
[]

[Kernels]
  # enforce c1 = (1-h(eta1))*c1l + h(eta1)*c1s
  [./PhaseConc1_eta1]
    type = KKSMultiPhaseConcentration
    variable = c1s
    cj = 'c1l c1s'
    hj_names = 'h1 h2'
    etas = 'eta1 eta2'
    c = c1
  [../]
  
  [./PhaseConc2_eta1]
    type = KKSMultiPhaseConcentration
    variable = c2s
    cj = 'c2l c2s'
    hj_names = 'h1 h2'
    etas = 'eta1 eta2'
    c = c2
  [../]
  
  [./PhaseConc1_eta2]
    type = KKSMultiPhaseConcentration
    variable = c1l
    cj = 'c1l c1s'
    hj_names = 'h1 h2'
    etas = 'eta1 eta2'
    c = c1
  [../]  
  [./PhaseConc2_eta2]
    type = KKSMultiPhaseConcentration
    variable = c2s
    cj = 'c2l c2s'
    hj_names = 'h1 h2'
    etas = 'eta1 eta2'
    c = c2
  [../] 


  # enforce pointwise equality of chemical potentials
  [./ChemPotSolute1]
    type = KKSPhaseChemicalPotential
    variable = c1l
    cb       = c1s
    fa_name  = fl
    fb_name  = fs
    args_a   = 'c2l'
    args_b   = 'c2s'
  [../]
  [./ChemPotSolute2]
    type = KKSPhaseChemicalPotential
    variable = c2l
    cb       = c2s
    fa_name  = fl
    fb_name  = fs
    args_a   = 'c1l'
    args_b   = 'c1s'

  [../]

  #
  # Cahn-Hilliard Equations
  #
  [./CHBulk1]
    type = KKSSplitCHCRes
    variable = c1
    ca       = c1l
    fa_name  = fl
    w        = w1
    args_a   = 'c2l'
  [../]
  [./CHBulk2]
    type = KKSSplitCHCRes
    variable = c2
    ca       = c2l
    fa_name  = fl
    w        = w2
    args_a   = 'c1l'
  [../]

  [./dc1dt]
    type = CoupledTimeDerivative
    variable = w1
    v = c1
  [../]
  [./dc2dt]
    type = CoupledTimeDerivative
    variable = w2
    v = c2
  [../]

  [./w1kernel]
    type = SplitCHWRes
    mob_name = M
    variable = w1
  [../]
  [./w2kernel]
    type = SplitCHWRes
    mob_name = M
    variable = w2
  [../]


  #
  # Allen-Cahn Equation
  #
  [./ACBulkF_eta1]
    type = KKSMultiACBulkF
    variable  = eta1
    Fj_names  = 'fl fs'
    hj_names  = 'h1 h2'
    gi_name   = g1
    eta_i     = eta1
    wi        = 1
    args      = 'c1l c1s c2l c2s eta2'
  [../]
  [./ACBulkC1_eta1]
      type = KKSMultiACBulkC
      variable  = eta1
      Fj_names  = 'fl fs'
      hj_names  = 'h1 h2'
      cj_names  = 'c1l c1s'
      eta_i     = eta1
      coupled_variables  = 'c2l c2s eta2'
    [../]  
    [./ACBulkC2_eta1]
      type = KKSMultiACBulkC
      variable  = eta1
      Fj_names  = 'fl fs'
      hj_names  = 'h1 h2'
      cj_names  = 'c2l c2s'
      eta_i     = eta1
      coupled_variables  = 'c1l c1s eta2'
    [../]  

  [./ACInterface_eta1]
    type = ACInterface
    variable = eta1
    kappa_name = eps_sq
  [../]
  [./detadt_eta1]
    type = TimeDerivative
    variable = eta1
  [../]
    [./ACBulkF_eta2]
      type = KKSMultiACBulkF
      variable  = eta2
      Fj_names  = 'fl fs'
      hj_names  = 'h1 h2'
      gi_name   = g2
      eta_i     = eta2
      wi        = 1
      args      = 'c1l c1s c2l c2s eta1'
    [../]
    [./ACBulkC1_eta2]
      type = KKSMultiACBulkC
      variable  = eta2
      Fj_names  = 'fl fs'
      hj_names  = 'h1 h2'
      cj_names  = 'c1l c1s'
      eta_i     = eta2
      coupled_variables  = 'c2l c2s eta1'
    [../]   
    [./ACBulkC2_eta2]
      type = KKSMultiACBulkC
      variable  = eta2
      Fj_names  = 'fl fs'
      hj_names  = 'h1 h2'
      cj_names  = 'c2l c2s'
      eta_i     = eta2
      coupled_variables  = 'c1l c1s eta1'
    [../]   
  
  [./ACInterface_eta2]
    type = ACInterface
    variable = eta2
    kappa_name = eps_sq
  [../]
  [./detadt_eta2]
    type = TimeDerivative
    variable = eta2
  [../]
[]

#[AuxKernels]
#  [./GlobalFreeEnergy]
#    variable = Fglobal
#    type = KKSGlobalFreeEnergy
#    fa_name = fl
#    fb_name = fs
#    w = 1.0
#  [../]
#[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      ilu          nonzero'

  l_max_its = 100
  nl_max_its = 100

  num_steps = 50
  dt = 1e-4
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

[Outputs]
  exodus = true
[]
