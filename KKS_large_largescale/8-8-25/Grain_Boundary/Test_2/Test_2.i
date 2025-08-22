[Mesh]
  [phasem_gr0]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 175
    ny = 175
    xmin = 0
    xmax = 350
    ymin = 0
    ymax = 350
    elem_type = QUAD4
  []
  [phasem_gr0_id]
    type = SubdomainIDGenerator
    input = phasem_gr0
    subdomain_id = 0
  []
  [phasem_gr1]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 175
    ny = 175
    xmin = 350
    xmax = 700
    ymin = 0
    ymax = 350
    elem_type = QUAD4
  []
  [phasem_gr1_id]
    type = SubdomainIDGenerator
    input = phasem_gr1
    subdomain_id = 1
  []
  [sticher]
    type = StitchedMeshGenerator
    inputs = 'phasem_gr0_id phasem_gr1_id'
    stitch_boundaries_pairs = 'right left'
    prevent_boundary_ids_overlap = false
  []
  [gb_sideset]
    type = SideSetsBetweenSubdomainsGenerator
    input = sticher
    primary_block = '0'
    paired_block  = '1'
    new_boundary  = 'gb'
  []
[]
[GlobalParams]
  op_num = 2
  var_name_base = gr
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Variables]
  [PolycrystalVariables]
  []
  # concentration
  [./c1]
    order = FIRST
    family = LAGRANGE
  [../]
  [./c2]
    order = FIRST
    family = LAGRANGE
  [../]

# phase concentration c1 in matrix
  [./c1m]
    order = FIRST
    family = LAGRANGE
  [../]
# phase concentration c1 in pv1
  [./c1pv1]
    order = FIRST
    family = LAGRANGE
  [../]
  # phase concentration c1 in pv2
  [./c1pv2]
    order = FIRST
    family = LAGRANGE
  [../]
 # phase concentration c1 in pv3
 [./c1pv3]
  order = FIRST
  family = LAGRANGE
[../]

# phase concentration c2 in matrix
  [./c2m]
    order = FIRST
    family = LAGRANGE
  [../]
  # phase concentration c2 in pv1
  [./c2pv1]
    order = FIRST
    family = LAGRANGE
  [../]
  # phase concentration c2 in pv2
  [./c2pv2]
    order = FIRST
    family = LAGRANGE
  [../]
 # phase concentration c2 in pv3
 [./c2pv3]
  order = FIRST
  family = LAGRANGE
[../]

  # order parameter m
  [./eta_m_gr0]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
  # order parameter pv1
  [./eta_pv1_gr0]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
  # order parameter pv2
  [./eta_pv2_gr0]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
# order parameter pv3
  [./eta_pv3_gr0]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
   # order parameter m
  [./eta_m_gr1]
    order = FIRST
    family = LAGRANGE
    block = 1
  [../]
  # order parameter pv1
  [./eta_pv1_gr1]
    order = FIRST
    family = LAGRANGE
    block = 1
  [../]
  # order parameter pv2
  [./eta_pv2_gr1]
    order = FIRST
    family = LAGRANGE
    block = 1
  [../]
# order parameter pv3
  [./eta_pv3_gr1]
    order = FIRST
    family = LAGRANGE
    block = 1
  [../]

[]

[ICs]
  [PolycrystalICs]
    [BicrystalBoundingBoxIC]
      x1 = 0
      y1 = 0
      x2 = 350
      y2 = 350
      block = '0 1'
    []
  []
  [./eta_pv1_gr0]
    variable = eta_pv1_gr0
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 324
    block = 0
  [../]
  [./eta_pv2_gr0]
    variable = eta_pv2_gr0
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 230	
    block = 0
  [../]
  [./eta_pv3_gr0]
    variable = eta_pv3_gr0
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 307	
    block = 0
  [../]
   [./eta_pv1_gr1]
    variable = eta_pv1_gr1
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 324
    block = 1
  [../]
  [./eta_pv2_gr1]
    variable = eta_pv2_gr1
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 230
    block = 1	
  [../]
  [./eta_pv3_gr1]
    variable = eta_pv3_gr1
    type = RandomIC
    min = -0.6
    max = 0.6
    seed = 307
    block = 1	
  [../]
  [./c1]
    variable = c1
    type = RandomIC
    min = 0.01775	
    max = 0.03025
    seed = 89	
  [../]
  [./c2]
    variable = c2
    type = RandomIC
    min = 0.03175	
    max = 0.04425
    seed = 89	
  [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution
    T = 500 # K
    wGB = 60 # nm
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 #Migration energy in eV
    GBenergy = 0.708 #GB energy in J/m^2
  [../]
[]

[Kernels]
 [./PolycrystalKernel]
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu            mumps            vinewtonrsls'

  l_max_its = 50
  nl_max_its = 25
  l_tol = 1.0e-3
  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-8

  end_time = 500 

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
  file_base = mesh
[]
