[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 2
  ny = 2
  nz = 2
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
    displacements = 'disp_x disp_y disp_z'
    component = 0
  [../]
  [./mech_y]
    type = LynxSolidMomentum
    variable = disp_y
    displacements = 'disp_x disp_y disp_z'
    component = 1
  [../]
  [./mech_z]
    type = LynxSolidMomentum
    variable = disp_z
    displacements = 'disp_x disp_y disp_z'
    component = 2
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
  [./vol_plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./yield]
    order = CONSTANT
    family = MONOMIAL
  [../]
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
  [./vol_plastic_strain_aux]
    type = LynxVolStrainAux
    variable = vol_plastic_strain
    strain_type = plastic
    execute_on = 'timestep_end'
  [../]
  [./yield_aux]
    type = MaterialRealAux
    variable = yield
    property = plastic_yield_function
    execute_on = 'timestep_end'
  [../]
[]

[BCs]
  [./no_x_left]
    type = PresetBC
    variable = disp_x
    boundary = left
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
  [./vx_comp]
    type = LynxVelocityBC
    variable = disp_x
    boundary = 'right'
    value = -1.6e-04
  [../]
  [./vy_comp]
    type = LynxVelocityBC
    variable = disp_y
    boundary = 'top'
    value = -1.6e-04
  [../]
  [./vz_comp]
    type = LynxVelocityBC
    variable = disp_z
    boundary = 'front'
    value = -1.6e-04
  [../]
[]

[Materials]
  [./deformation]
    type = LynxDamageDeformation
    displacements = 'disp_x disp_y disp_z'
    bulk_modulus = 2.0e+09
    shear_modulus = 2.0e+09
    friction_angle = 30
    critical_pressure = 8.0e+07
    plastic_viscosity = 5.0e+10
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
  [./Ev_p]
    type = ElementAverageValue
    variable = vol_plastic_strain
    outputs = csv
  [../]
[]

[Preconditioning]
  [./precond]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -ksp_max_it -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres asm 1E-05 1E-10 200 500 lu NONZERO'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  start_time = 0.0
  end_time = 250
  dt = 2.5
[]

[Outputs]
  execute_on = 'timestep_end'
  print_linear_residuals = true
  perf_graph = true
  exodus = true
  csv = true
[]

