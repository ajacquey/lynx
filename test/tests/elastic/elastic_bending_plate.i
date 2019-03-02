[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 35
  ny = 4
  nz = 4
  xmin = 0
  xmax = 10e+03
  ymin = -0.5e+03
  ymax = 0.5e+03
  zmin = -0.5e+03
  zmax = 0.5e+03
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
  [./Pe]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ev]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./S_dev]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Sxz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Exz]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./pe_aux]
    type = LynxEffectivePressureAux
    variable = Pe
  [../]
  [./Ev_aux]
    type = LynxVolStrainAux
    variable = Ev
  [../]
  [./S_aux]
    type = LynxVonMisesStressAux
    variable = S_dev
  [../]
  [./Sxz_aux]
    type = LynxStressAux
    variable = Sxz
    index_i = 0
    index_j = 2
  [../]
  [./Exz_aux]
    type = LynxStrainAux
    variable = Exz
    index_i = 0
    index_j = 2
  [../]
[]

[BCs]
  [./no_ux]
    type = PresetBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  [../]
  [./no_uy]
    type = PresetBC
    variable = disp_y
    boundary = 'left right bottom top front back'
    value = 0.0
  [../]
  [./no_uz]
    type = PresetBC
    variable = disp_z
    boundary = 'left'
    value = 0.0
  [../]
[]

[Materials]
  [./elastic_mat]
    type = LynxDeformation
    bulk_modulus = 6.6e+10
    shear_modulus = 3.95e+10
  [../]
  [./density]
    type = LynxDensityConstant
    solid_density = 100
    has_gravity = true
  [../]
[]

[Preconditioning]
  [./precond]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -ksp_max_it -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres asm 1E-10 1E-10 200 500 lu NONZERO'
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
  #execute_on = 'timestep_end'
  print_linear_residuals = true
  perf_graph = true
  exodus = true
[]