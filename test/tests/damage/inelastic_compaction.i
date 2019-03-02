[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 4
  nz = 4
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  zmin = 0
  zmax = 1
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./mech_x]
    type = LynxSolidMomentum
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./mech_y]
    type = LynxSolidMomentum
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./mech_z]
    type = LynxSolidMomentum
    variable = disp_z
    component = 2
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[AuxVariables]
  [./pressure]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vol_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vol_inelastic_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vol_inelastic_strain_rate]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_ratio]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./damage_yield]
    order = CONSTANT
    family = MONOMIAL
  [../]
#  [./porosity_yield]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
[]

[AuxKernels]
  [./pressure_aux]
    type = LynxEffectivePressureAux
    variable = pressure
    execute_on = 'timestep_end'
  [../]
  [./vol_strain_aux]
    type = LynxVolStrainAux
    variable = vol_strain
    execute_on = 'timestep_end'
  [../]
  [./vol_inelastic_strain_aux]
    type = LynxVolStrainAux
    variable = vol_inelastic_strain
    strain_type = inelastic
    execute_on = 'timestep_end'
  [../]
  [./vol_inelastic_strain_rate_aux]
    type = LynxVolStrainRateAux
    variable = vol_inelastic_strain_rate
    strain_type = inelastic
    execute_on = 'timestep_end'
  [../]
  [./strain_ratio_aux]
    type = LynxStrainRatioAux
    variable = strain_ratio
    execute_on = 'timestep_end'
  [../]
  [./damage_yield_aux]
    type = MaterialRealAux
    variable = damage_yield
    property = damage_yield_function
  [../]
#  [./porosity_yield_aux]
#    type = MaterialRealAux
#    variable = porosity_yield
#    property = porosity_yield_function
#  [../]
[]

[BCs]
  [./no_x_left]
    type = PresetBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  [../]
  [./no_y_bottom]
    type = PresetBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  [../]
  [./no_z_back]
    type = PresetBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  [../]
  [./compression_x]
    type = LynxVelocityBC
    variable = disp_x
    boundary = 'right'
    value = -1.0e-12
  [../]
  [./compression_y]
    type = LynxVelocityBC
    variable = disp_y
    boundary = 'top'
    value = -1.0e-12
  [../]
  [./compression_z]
    type = LynxVelocityBC
    variable = disp_z
    boundary = 'front'
    value = -1.0e-12
  [../]
[]

[Materials]
  [./deformation]
    type = LynxDamageDeformation
    displacements = 'disp_x disp_y disp_z'
    bulk_modulus = 1.0e+10
    shear_modulus = 1.0e+10
    damage_modulus = 0.0
    friction_angle = 30
    porous_coupling = 10.0
    plastic_viscosity = 5.0e+20
    [../]
[]

[Postprocessors]
  [./P]
    type = ElementAverageValue
    variable = pressure
    outputs = csv
  [../]
  [./Ev]
    type = ElementAverageValue
    variable = vol_strain
    outputs = csv
  [../]
  [./Ev_in]
    type = ElementAverageValue
    variable = vol_inelastic_strain
    outputs = csv
  [../]
  [./E_dot_v_in]
    type = ElementAverageValue
    variable = vol_inelastic_strain_rate
    outputs = csv
  [../]
  [./Xi]
    type = ElementAverageValue
    variable = strain_ratio
    outputs = csv
  [../]
  [./F]
    type = ElementAverageValue
    variable = damage_yield
    outputs = csv
  [../]
#  [./Fp]
#    type = ElementAverageValue
#    variable = porosity_yield
#    outputs = csv
#  [../]
[]

[Preconditioning]
  [./precond]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -ksp_max_it -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres asm 1E-00 1E-15 200 500 ilu NONZERO'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0.0
  end_time = 3.1536e+11
  dt = 6.3072e+09
[]

[Outputs]
  execute_on = 'timestep_end'
  print_linear_residuals = false
  perf_graph = true
  exodus = true
  csv = true
[]
