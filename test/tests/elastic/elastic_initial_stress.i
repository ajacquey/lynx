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

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[Kernels]
  [./mech_x]
    type = LynxSolidMomentum
    variable = disp_x
    component = 0
  [../]
  [./mech_y]
    type = LynxSolidMomentum
    variable = disp_y
    component = 1
  [../]
  [./mech_z]
    type = LynxSolidMomentum
    variable = disp_z
    component = 2
  [../]
[]

[AuxVariables]
  [./Sxx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Syy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Szz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Exx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Eyy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ezz]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./Sxx_aux]
    type = LynxStressAux
    variable = Sxx
    index_i = 0
    index_j = 0
    execute_on = 'timestep_end'
  [../]
  [./Syy_aux]
    type = LynxStressAux
    variable = Syy
    index_i = 1
    index_j = 1
    execute_on = 'timestep_end'
  [../]
  [./Szz_aux]
    type = LynxStressAux
    variable = Szz
    index_i = 2
    index_j = 2
    execute_on = 'timestep_end'
  [../]
  [./Exx_aux]
    type = LynxStrainAux
    variable = Exx
    index_i = 0
    index_j = 0
    execute_on = 'timestep_end'
  [../]
  [./Eyy_aux]
    type = LynxStrainAux
    variable = Eyy
    index_i = 1
    index_j = 1
    execute_on = 'timestep_end'
  [../]
  [./Ezz_aux]
    type = LynxStrainAux
    variable = Ezz
    index_i = 2
    index_j = 2
    execute_on = 'timestep_end'
  [../]
[]

[BCs]
  [./no_ux]
    type = PresetBC
    variable = disp_x
    boundary = 'left right'
    value = 0.0
  [../]
  [./no_uy]
    type = PresetBC
    variable = disp_y
    boundary = 'bottom top'
    value = 0.0
  [../]
  [./no_uz]
    type = PresetBC
    variable = disp_z
    boundary = 'back front'
    value = 0.0
  [../]
[]

[Materials]
  [./elastic_mat]
    type = LynxDeformation
    bulk_modulus = 1.0e+09
    shear_modulus = 1.0e+09
    initial_stress = '-40.0e+06 0.0 0.0 -50.0e+06 0.0 -100.0e+06'
  [../]
[]

[Preconditioning]
  [./precond]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -ksp_max_it -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres asm 1E-05 1E-10 20 50 ilu NONZERO'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  start_time = 0.0
  end_time = 1.0
  dt = 1.0
[]

[Outputs]
  execute_on = 'timestep_end'
  print_linear_residuals = true
  print_perf_log = true
  exodus = true
[]