# This test is for the 3-phase KKS model


[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 150
  ny = 150
#  nz = 2
  xmin = 0
  xmax = 300
  ymin = 0
  ymax = 300
  zmin = 0
  zmax = 0
#  elem_type = QUAD4
[]

[BCs]
  [./u_right_pull]
    type = Pressure
    displacements = 'disp_x disp_y'
    variable = disp_x
    boundary = right
    factor = 0
  [../]
  [./all_Al]
    type =  NeumannBC
    variable = 'Al'
    boundary = 'left right top bottom'
    value = 0
  [../]	
  [./all_Nb]
    type =  NeumannBC
    variable = 'Nb'
    boundary = 'left right top bottom'
    value = 0
  [../]	
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
[]


[AuxVariables]
  [./Energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./e_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Bounds]
#  [./eta_m_upper_bound]
#    type = ConstantBoundsAux
#    variable = bounds_dummy
#    bounded_variable = eta_m
#    bound_type = upper
#    bound_value = 1
#  [../]
#  [./eta_m_lower_bound]
#    type = ConstantBoundsAux
#    variable = bounds_dummy
#    bounded_variable = eta_m
#    bound_type = lower
#    bound_value = 0
#  [../]

  [./eta_pv1_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv1
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv1_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv1
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_pv2_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv2
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv2_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv2
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_pv3_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv3
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv3_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv3
    bound_type = lower
    bound_value = -1
  [../]
  [./eta_pv4_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv4
    bound_type = upper
    bound_value = 1
  [../]
  [./eta_pv4_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = eta_pv4
    bound_type = lower
    bound_value = -1
  [../]
[]


[Variables]
  # concentration
  [./Al]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Nb]
    order = FIRST
    family = LAGRANGE
  [../]

# phase concentration Al in matrix
  [./Alm]
    order = FIRST
    family = LAGRANGE
  [../]
# phase concentration Al in pv1
  [./Alpv1]
    order = FIRST
    family = LAGRANGE
  [../]
  # phase concentration Al in pv2
  [./Alpv2]
    order = FIRST
    family = LAGRANGE
  [../]
 # phase concentration Al in pv3
 [./Alpv3]
  order = FIRST
  family = LAGRANGE
[../]
# phase concentration Al in pv4
[./Alpv4]
  order = FIRST
  family = LAGRANGE
[../]
# phase concentration Nb in matrix
  [./Nbm]
    order = FIRST
    family = LAGRANGE
  [../]
  # phase concentration Nb in pv1
  [./Nbpv1]
    order = FIRST
    family = LAGRANGE
  [../]
  # phase concentration Nb in pv2
  [./Nbpv2]
    order = FIRST
    family = LAGRANGE
  [../]
 # phase concentration Nb in pv3
 [./Nbpv3]
  order = FIRST
  family = LAGRANGE
[../]
# phase concentration Nb in pv4
[./Nbpv4]
  order = FIRST
  family = LAGRANGE
[../]

# order parameter m
[./eta_m]
  order = FIRST
  family = LAGRANGE
[../]
# order parameter pv1
[./eta_pv1]
  order = FIRST
  family = LAGRANGE
[../]
  # order parameter pv2
[./eta_pv2]
  order = FIRST
  family = LAGRANGE
[../]
# order parameter pv3
[./eta_pv3]
  order = FIRST
  family = LAGRANGE
[../]
# order parameter pv4
[./eta_pv4]
  order = FIRST
  family = LAGRANGE
[../]

  [./disp_x]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
[]

[Functions]
  [./bc_func]
    type = ParsedFunction
    value = sin(alpha*pi*x)
    vars = alpha
    vals = 16
  [../]
[]

[ICs]
  [./eta_pv1]
    variable = eta_pv1
    type = RandomIC
    min = -0.1625
    max = 0.1625
    seed = 192
  [../]
  [./eta_pv2]
    variable = eta_pv2
    type = RandomIC
    min = -0.1625
    max = 0.1625
    seed = 192	
  [../]
  [./eta_pv3]
    variable = eta_pv3
    type = RandomIC
    min = -0.1625
    max = 0.1625
    seed = 192	
  [../]
  [./eta_pv4]
    variable = eta_pv4
    type = RandomIC
    min = -0.1625
    max = 0.1625
    seed = 192	
  [../]
  [./Al]
    variable = Al
    type = RandomIC
    min = 0.01775	
    max = 0.03025
    seed = 89	
  [../]
  [./Nb]
    variable = Nb
    type = RandomIC
    min = 0.03175	
    max = 0.04425
    seed = 89	
  [../]

[]


[Materials]
  [./fm]
    type = DerivativeParsedMaterial
    f_name = fc_m
    args = 'Alm Nbm'
    function = '50.0*((Alm-0.0161)^2+2*(Nbm-0.00723)^2)'
  [../]
  # Elastic energy of the phase 0
  [./elastic_free_energy_m]
    type = ElasticEnergyMaterial
    base_name = phasem
    f_name = fe_m
    args = ' '
  [../]
  # Total free energy of the phase 0
  [./Total_energy_m]
    type = DerivativeSumMaterial
    f_name = Fm
    sum_materials = 'fc_m fe_m'
    args = 'Alm Nbm'
  [../]

  #gamma_prime
  [./fc_pv1]
    type = DerivativeParsedMaterial
    f_name = fc_pv1
    args = 'Alpv1 Nbpv1'
    function =  '50.0*((Alpv1-0.187)^2+2*(Nbpv1-0.0157)^2)'
  [../]
  # Elastic energy of the phase 1
  [./elastic_free_energy_pv1]
    type = ElasticEnergyMaterial
    base_name = phasepv1
    f_name = fe_pv1
    args = ' '
  [../]
  # Total free energy of the phase 1
  [./Total_energy_pv1]
    type = DerivativeSumMaterial
    f_name = Fpv1
    sum_materials = 'fc_pv1 fe_pv1'
    args = 'Alpv1 Nbpv1'
  [../]

  #gamma_double prime variant1
  [./f2]
    type = DerivativeParsedMaterial
    f_name = fc_pv2
    args = 'Alpv2 Nbpv2'
    function = '50.0*((Alpv2-0.00727)^2+2*(Nbpv2-0.196)^2)'
  [../]
  # Elastic energy of the phase 2
  [./elastic_free_energy_pv2]
    type = ElasticEnergyMaterial
    base_name = phasepv2
    f_name = fe_pv2
    args = ' '
  [../]
  # Total free energy of the phase 2
  [./Total_energy_pv2]
    type = DerivativeSumMaterial
    f_name = Fpv2
    sum_materials = 'fc_pv2 fe_pv2'
    args = 'Alpv2 Nbpv2'
  [../]

#gamma_double prime variant2
  [./f3]
    type = DerivativeParsedMaterial
    f_name = fc_pv3
    args = 'Alpv3 Nbpv3'
    function = '50.0*((Alpv3-0.00727)^2+2*(Nbpv3-0.196)^2)'
  [../]
    # Elastic energy of the phase 3
  [./elastic_free_energy_pv3]
    type = ElasticEnergyMaterial
    base_name = phasepv3
    f_name = fe_pv3
    args = ' '
  [../]
    # Total free energy of the phase 3
  [./Total_energy_pv3]
    type = DerivativeSumMaterial
    f_name = Fpv3
    sum_materials = 'fc_pv3 fe_pv3'
    args = 'Alpv3 Nbpv3'
  [../]
#gamma_double prime variant3
  [./f4]
    type = DerivativeParsedMaterial
    f_name = fc_pv4
    args = 'Alpv4 Nbpv4'
    function = '50.0*((Alpv4-0.00727)^2+2*(Nbpv4-0.196)^2)'
  [../]
  # Elastic energy of the phase 4
  [./elastic_free_energy_pv4]
    type = ElasticEnergyMaterial
    base_name = phasepv4
    f_name = fe_pv4
    args = ' '
  [../]
    # Total free energy of the phase 4
  [./Total_energy_pv4]
    type = DerivativeSumMaterial
    f_name = Fpv4
    sum_materials = 'fc_pv4 fe_pv4'
    args = 'Alpv4 Nbpv4'
  [../]

  # Switching functions for each phase
  # hm(eta_pv1, eta_pv2, eta_m)
  [./hm]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_m
    all_etas = 'eta_pv1 eta_pv2 eta_pv3 eta_m'
    h_name = hm
  [../]
  # hpv1(eta_pv1 eta_pv2 eta_pv3 eta_pv4 eta_m)
  [./hpv1]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv1
    all_etas = 'eta_pv1 eta_pv2 eta_pv3 eta_pv4 eta_m'
    h_name = hpv1
  [../]
  # hpv2(eta_pv1 eta_pv2 eta_pv3 eta_pv4 eta_m)
  [./hpv2]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv2
    all_etas = 'eta_pv1 eta_pv2 eta_pv3 eta_pv4 eta_m'
    h_name = hpv2
  [../]
  # hpv3(eta_pv1 eta_pv2 eta_pv3 eta_pv4 eta_m)
[./hpv3]
    type = SwitchingFunctionMultiPhaseMaterial
    phase_etas = eta_pv3
    all_etas = 'eta_pv1 eta_pv2 eta_pv3 eta_pv4 eta_m'
    h_name = hpv3
[../]
# hpv4(eta_pv1 eta_pv2 eta_pv3 eta_pv4 eta_m)
[./hpv4]
  type = SwitchingFunctionMultiPhaseMaterial
  phase_etas = eta_pv4
  all_etas = 'eta_pv1 eta_pv2 eta_pv3 eta_pv4 eta_m'
  h_name = hpv4
[../]
  # Coefficients for diffusion equation
  [./Dhm]
    type = DerivativeParsedMaterial
    material_property_names = 'D hm'
    function = D*hm
    f_name = Dhm
  [../]
  [./Dhpv1]
    type = DerivativeParsedMaterial
    material_property_names = 'D hpv1'
    function = D*hpv1
    f_name = Dhpv1
  [../]
  [./Dhpv2]
    type = DerivativeParsedMaterial
    material_property_names = 'D hpv2'
    function = D*hpv2
    f_name = Dhpv2
  [../]
  [./Dhpv3]
    type = DerivativeParsedMaterial
    material_property_names = 'D hpv3'
    function = D*hpv3
    f_name = Dhpv3
  [../]
  [./Dhpv4]
    type = DerivativeParsedMaterial
    material_property_names = 'D hpv4'
    function = D*hpv4
    f_name = Dhpv4
  [../]

# Barrier functions for each phase
  [./gm]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta_m
    function_name = gm
  [../]
  [./gpv1]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta_pv1
    function_name = gpv1
  [../]
  [./gpv2]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta_pv2
    function_name = gpv2
  [../]
  [./gpv3]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta_pv3
    function_name = gpv3
  [../]
  [./gpv4]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta_pv4
    function_name = gpv4
  [../]

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'L    kappa  D  misfit   W'
    prop_values = '0.3  0.01   1  0.005   0.01'
  [../]

  #Mechanical properties
  [./Stiffness_phasem]
    type = ComputeElasticityTensor
    C_ijkl = '2008.5 448.9 918.5 2008.5 918.5 1538.7 779.8 779.8 309.9'    
    base_name = phasem
    fill_method = symmetric9
  [../]
  [./Stiffness_phasepv1]
    type = ComputeElasticityTensor
    C_ijkl = '1982 534 496 1606 605 1788 985 1056 388'    
    base_name = phasepv1
    fill_method = symmetric9
  [../]
  [./Stiffness_phasepv2]
    type = ComputeElasticityTensor
    C_ijkl = '1606 534 605 1982 496 1788 1056 985 388'
    base_name = phasepv2
    fill_method = symmetric9
  [../]
  [./Stiffness_phasepv3]
    type = ComputeElasticityTensor
    C_ijkl = '1982 534 496 1606 605 1788 985 1056 388'
    base_name = phasepv3
    fill_method = symmetric9
  [../]
  [./Stiffness_phasepv4]
    type = ComputeElasticityTensor
    C_ijkl = '1982 534 496 1606 605 1788 985 1056 388'
      base_name = phasepv4
      fill_method = symmetric9
  [../]

  [./stress_phasepv1]
    type = ComputeLinearElasticStress
    base_name = phasepv1
  [../]
  [./stress_phasepv2]
    type = ComputeLinearElasticStress
    base_name = phasepv2
  [../]
  [./stress_phasepv3]
    type = ComputeLinearElasticStress
    base_name = phasepv3
  [../]
  [./stress_phasepv4]
    type = ComputeLinearElasticStress
    base_name = phasepv4
  [../]
  [./stress_phasem]
    type = ComputeLinearElasticStress
    base_name = phasem
  [../]

  [./strain_phasem]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasem
  [../]
  [./strain_phasepv1]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv1
    eigenstrain_names = eigenstrainpv1
  [../]
  [./strain_phasepv2]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv2
    eigenstrain_names = eigenstrainpv2
  [../]
  [./strain_phasepv3]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv3
    eigenstrain_names = eigenstrainpv3
  [../]
  [./strain_phasepv4]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    base_name = phasepv4
    eigenstrain_names = eigenstrainpv4
  [../]

  [./eigen_strainpv1]
    type = ComputeEigenstrain
    base_name = phasepv1
    eigen_base = '1 0 0 0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv1
  [../]

  [./eigen_strainpv2]
    type = ComputeEigenstrain
    base_name = phasepv2
    eigen_base = '0 1 0 0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv2
  [../]
  [./eigen_strainpv3]
    type = ComputeEigenstrain
    base_name = phasepv3
    eigen_base = '1 0 0 0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv3
  [../]
  [./eigen_strainpv4]
    type = ComputeEigenstrain
    base_name = phasepv4
    eigen_base = '1 0 0 0 0 0'
    prefactor = misfit
    eigenstrain_name = eigenstrainpv4
  [../]


  # Generate the global stress from the phase stresses
  [./global_stress]
    type = MultiPhaseStressMaterial
    phase_base = 'phasepv1 phasepv2 phasepv3 phasepv4 phasem'
    h          = 'hpv1     hpv2   hpv3 hpv4  hm'
  [../]

  [./global_strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
  
  # Kernels for Allen-Cahn equation for eta_pv1
  [./deta_pv1_dt]
    type = TimeDerivative
    variable = eta_pv1
  [../]
  [./ACBulkFpv1]
    type = KKSMultiACBulkF
    variable  = eta_pv1
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    gi_name   = gpv1
    eta_i     = eta_pv1
    wi        = 0.01
    args      = 'Alpv1 Alpv2 Alpv3 Alpv4 Alm Nbpv1 Nbpv2 Nbpv3 Nbpv4 Nbm eta_pv1 eta_pv3 eta_pv4 eta_m'
  [../]
  [./ACBulkCpv1_Al]
    type = KKSMultiACBulkC
    variable  = eta_pv1
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'Alpv1 Alpv2 Alpv3 Alpv4 Alm'
    eta_i     = eta_pv1
    args      = 'eta_pv2 eta_pv3 eta_pv4 eta_m'
  [../]
  [./ACBulkCpv1_Nb]
    type = KKSMultiACBulkC
    variable  = eta_pv1
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'Nbpv1 Nbpv2 Nbpv3 Nbpv4 Nbm'
    eta_i     = eta_pv1
    args      = 'eta_pv2 eta_pv3 eta_pv4 eta_m'
  [../]
  [./ACInterfacepv1]
    type = ACInterface
    variable = eta_pv1
    kappa_name = kappa
  [../]

  # Kernels for Allen-Cahn equation for eta_pv2
  [./deta_pv2_dt]
    type = TimeDerivative
    variable = eta_pv2
  [../]
  [./ACBulkFpv2]
    type = KKSMultiACBulkF
    variable  = eta_pv2
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    gi_name   = gpv2
    eta_i     = eta_pv2
    wi        = 0.01
    args      = 'Alpv1 Alpv2 Alpv3 Alpv4 Alm Nbpv1 Nbpv2 Nbpv3 Nbpv4 Nbm eta_pv1 eta_pv3 eta_pv4 eta_m'
  [../]
  [./ACBulkCpv2_Al]
    type = KKSMultiACBulkC
    variable  = eta_pv2
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'Alpv1 Alpv2 Alpv3 Alpv4 Alm'
    eta_i     = eta_pv2
    args      = 'eta_pv1 eta_pv3 eta_pv4 eta_m'
  [../]
  [./ACBulkCpv2_Nb]
    type = KKSMultiACBulkC
    variable  = eta_pv2
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'Nbpv1 Nbpv2 Nbpv3 Nbpv4 Nbm'
    eta_i     = eta_pv2
    args      = 'eta_pv1 eta_pv3 eta_pv4 eta_m'
  [../]
  [./ACInterfacepv2]
    type = ACInterface
    variable = eta_pv2
    kappa_name = kappa
  [../]

  # Kernels for Allen-Cahn equation for eta_pv3
  [./deta_pv3_dt]
    type = TimeDerivative
    variable = eta_pv3
  [../]
  [./ACBulkFpv3]
    type = KKSMultiACBulkF
    variable  = eta_pv3
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    gi_name   = gpv3
    eta_i     = eta_pv3
    wi        = 0.01
    args      = 'Alpv1 Alpv2 Alpv3 Alpv4 Alm Nbpv1 Nbpv2 Nbpv3 Nbpv4 Nbm eta_pv1 eta_pv2 eta_pv4 eta_m'
  [../]
  [./ACBulkCpv3_Al]
    type = KKSMultiACBulkC
    variable  = eta_pv3
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'Alpv1 Alpv2 Alpv3 Alpv4 Alm'
    eta_i     = eta_pv3
    args      = 'eta_pv1 eta_pv2 eta_pv4 eta_m'
  [../]
  [./ACBulkCpv3_Nb]
    type = KKSMultiACBulkC
    variable  = eta_pv3
    Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
    cj_names  = 'Nbpv1 Nbpv2 Nbpv3 Nbpv4 Nbm'
    eta_i     = eta_pv3
    args      = 'eta_pv1 eta_pv2 eta_pv4 eta_m'
  [../]
  [./ACInterfacepv3]
    type = ACInterface
    variable = eta_pv3
    kappa_name = kappa
  [../]
 # Kernels for Allen-Cahn equation for eta_pv4
 [./deta_pv4_dt]
  type = TimeDerivative
  variable = eta_pv4
[../]
[./ACBulkFpv4]
  type = KKSMultiACBulkF
  variable  = eta_pv4
  Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
  hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
  gi_name   = gpv4
  eta_i     = eta_pv4
  wi        = 0.01
  args      = 'Alpv1 Alpv2 Alpv3 Alpv4 Alm Nbpv1 Nbpv2 Nbpv3 Nbpv4 Nbm eta_pv1 eta_pv2 eta_pv3 eta_m'
[../]
[./ACBulkCpv4_Al]
  type = KKSMultiACBulkC
  variable  = eta_pv4
  Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
  hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
  cj_names  = 'Alpv1 Alpv2 Alpv3 Alpv4 Alm'
  eta_i     = eta_pv4
  args      = 'eta_pv1 eta_pv2 eta_pv3 eta_m'
[../]
[./ACBulkCpv4_Nb]
  type = KKSMultiACBulkC
  variable  = eta_pv4
  Fj_names  = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
  hj_names  = 'hpv1 hpv2 hpv3 hpv4 hm'
  cj_names  = 'Nbpv1 Nbpv2 Nbpv3 Nbpv4 Nbm'
  eta_i     = eta_pv4
  args      = 'eta_pv1 eta_pv2 eta_pv3 eta_m'
[../]
[./ACInterfacepv4]
  type = ACInterface
  variable = eta_pv4
  kappa_name = kappa
[../]

# Kernels for constraint equation |eta_pv1| + |eta_pv2| + eta_m = 1
  # eta3 is the nonlinear variable for the constraint equation
  [./eta_mreaction]
    type = MatReaction
    variable = eta_m
    reaction_rate = 1
 [../]
  [./eta_pv1reaction]
    type = MatReaction
    variable = eta_m
    v = eta_pv1
    reaction_rate = 1
  [../]
  [./eta_pv2reaction]
    type = MatReaction
    variable = eta_m
    v = eta_pv2
    reaction_rate = 1
  [../]
  [./eta_pv3reaction]
    type = MatReaction
    variable = eta_m
    v = eta_pv3
    reaction_rate = 1
  [../]
  [./eta_pv4reaction]
    type = MatReaction
    variable = eta_m
    v = eta_pv4
    reaction_rate = 1
  [../]
  [./one]
    type = BodyForce
    variable = eta_m
    value = -1.0
  [../]

  #Kernels for diffusion equation of Al
  [./diff_time_Al]
    type = TimeDerivative
    variable = Al
  [../]
  [./diff_Alm]
    type = MatDiffusion
    variable = Al
    diffusivity = Dhm
    v = Alm
  [../]
  [./diff_Alpv1]
    type = MatDiffusion
    variable = Al
    diffusivity = Dhpv1
    v = Alpv1
  [../]
  [./diff_Alpv2]
    type = MatDiffusion
    variable = Al
    diffusivity = Dhpv2
    v = Alpv2
  [../]
  [./diff_Alpv3]
    type = MatDiffusion
    variable = Al
    diffusivity = Dhpv3
    v = Alpv3
  [../]
  [./diff_Alpv4]
    type = MatDiffusion
    variable = Al
    diffusivity = Dhpv4
    v = Alpv4
  [../]

  #Kernels for diffusion equation of Nb
  [./diff_time_Nb]
    type = TimeDerivative
    variable = Nb
  [../]
  [./diff_Nbm]
    type = MatDiffusion
    variable = Nb
    diffusivity = Dhm
    v = Nbm
  [../]
  [./diff_Nbpv1]
    type = MatDiffusion
    variable = Nb
    diffusivity = Dhpv1
    v = Nbpv1
  [../]
  [./diff_Nbpv2]
    type = MatDiffusion
    variable = Nb
    diffusivity = Dhpv2
    v = Nbpv2
  [../]
  [./diff_Nbpv3]
    type = MatDiffusion
    variable = Nb
    diffusivity = Dhpv3
    v = Nbpv3
  [../]
  [./diff_Nbpv4]
    type = MatDiffusion
    variable = Nb
    diffusivity = Dhpv4
    v = Nbpv4
  [../]    
    
  # Phase concentration constraints
  [./chempot1m_pv1]
    type = KKSPhaseChemicalPotential
    variable = Alm
    cb       = Alpv1
    fa_name  = Fm
    fb_name  = Fpv1
  [../]
  [./chempot1m_pv2]
    type = KKSPhaseChemicalPotential
    variable = Alpv2
    cb       = Alm
    fa_name  = Fpv2
    fb_name  = Fm
  [../]
  [./chempot1m_pv3]
    type = KKSPhaseChemicalPotential
    variable = Alpv3
    cb       = Alm
    fa_name  = Fpv3
    fb_name  = Fm
  [../]
  [./chempot1m_pv4]
    type = KKSPhaseChemicalPotential
    variable = Alpv4
    cb       = Alm
    fa_name  = Fpv4
    fb_name  = Fm
  [../]
  [./chempot2m_pv1]
    type = KKSPhaseChemicalPotential
    variable = Nbm
    cb       = Nbpv1
    fa_name  = Fm
    fb_name  = Fpv1
  [../]
  [./chempot2m_pv2]
    type = KKSPhaseChemicalPotential
    variable = Nbpv2
    cb       = Nbm
    fa_name  = Fpv2
    fb_name  = Fm
  [../]
  [./chempot2m_pv3]
    type = KKSPhaseChemicalPotential
    variable = Nbpv3
    cb       = Nbm
    fa_name  = Fpv3
    fb_name  = Fm
  [../]
  [./chempot2m_pv4]
    type = KKSPhaseChemicalPotential
    variable = Nbpv4
    cb       = Nbm
    fa_name  = Fpv4
    fb_name  = Fm
  [../]
    
  [./phaseconcentration_Alpv1]
    type = KKSMultiPhaseConcentration
    variable = Alpv1
    cj = 'Alm Alpv1 Alpv2 Alpv3 Alpv4'
    hj_names = 'hm hpv1 hpv2 hpv3 hpv4'
    etas = 'eta_m eta_pv1 eta_pv2 eta_pv3 eta_pv4'
    c = Al
  [../]
  [./phaseconcentration_Nbpv1]
    type = KKSMultiPhaseConcentration
    variable = Nbpv1
    cj = 'Nbm Nbpv1 Nbpv2 Nbpv3 Nbpv4'
    hj_names = 'hm hpv1 hpv2 hpv3 hpv4'
    etas = 'eta_m eta_pv1 eta_pv2 eta_pv3 eta_pv4'
    c = Nb
  [../]
  

[]

[AuxKernels]
  [./Energy_total]
    type = KKSMultiFreeEnergy
    Fj_names = 'Fpv1 Fpv2 Fpv3 Fpv4 Fm'
    hj_names = 'hpv1 hpv2 hpv3 hpv4 hm'
    gj_names = 'gpv1 gpv2 gpv3 gpv4 gm'
    variable = Energy
    w = 1
    interfacial_vars =  'eta_pv1  eta_pv2 eta_pv3 eta_pv4 eta_m'
    kappa_names =       'kappa kappa kappa kappa kappa'
  [../]
  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
    execute_on = timestep_end
    block = 0
  [../]
#  [./e_xx]
#    type = RankTwoAux
#    variable = e_xx
#    rank_two_tensor = phase1_lage
#    index_j = 0
#    index_i = 0
#    execute_on = timestep_end
#    block = 0
#  [../]

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

  end_time = 50000

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 5e-4
    cutback_factor = 0.75
    growth_factor = 1.2
    optimal_iterations = 20
  [../]

[]

[Preconditioning]
  active = 'full'
  [./full]
    type = SMP
    full = true
  [../]
  [./mydebug]
    type = FDP
    full = true
  [../]
[]

[Postprocessors]
   [./dofs]
     type = NumDOFs
   [../]
   [./h1_error]
     type = ElementH1Error
     variable = eta_pv1
     function = bc_func
   [../]
[]

[Outputs]
  exodus = true
  [./table]
    type = CSV
    execute_on = timestep_end
  [../]
  [./checkpoint]
     type = Checkpoint
     num_files = 10
     interval = 10
  [../]
[]
