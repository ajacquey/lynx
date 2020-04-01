[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 10
  xmin = 0
  xmax = 0.1
  ymin = 0
  ymax = 0.1
  zmin = 0
  zmax = 1
[]

[Variables]
  [./pf]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./darcy_div]
    type = LynxADHydroDarcy
    variable = pf
  [../]
[]

[AuxVariables]
  [./darcy_vel_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./darcy_vel_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./darcy_vel_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./darcy_vel_x_aux]
    type = LynxDarcyVelocityAux
    variable = darcy_vel_x
    component = 0
    fluid_pressure = pf
  [../]
  [./darcy_vel_y_aux]
    type = LynxDarcyVelocityAux
    variable = darcy_vel_y
    component = 1
    fluid_pressure = pf
  [../]
  [./darcy_vel_z_aux]
    type = LynxDarcyVelocityAux
    variable = darcy_vel_z
    component = 2
    fluid_pressure = pf
  [../]
[]

[BCs]
  [./top_p0]
    type = ADDirichletBC
    variable = pf
    boundary = 'front'
    value = 0.0
  [../]
  [./bottom_p1]
    type = ADDirichletBC
    variable = pf
    boundary = 'back'
    value = 2.0
  [../]
[]

[Materials]
  [./hydro_const]
    type = LynxADHydroConstant
    fluid_viscosity = 1.0
    permeability = 0.5
    porosity = 0.1
  [../]
  [./density_const]
    type = LynxADDensityConstant
    has_gravity = true
    gravity_acceleration = 1.0
    fluid_density = 1.0
    solid_density = 0.0
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
  type = Steady
  solve_type = 'NEWTON'
  automatic_scaling = true
[]

[Outputs]
  perf_graph = true
  exodus = true
[]
