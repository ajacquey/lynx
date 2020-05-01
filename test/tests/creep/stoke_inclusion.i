[Mesh]
  type = FileMesh
  file = stoke_inclusion.msh
  boundary_id = '1 2 3 4'
  boundary_name = 'bottom right top left'
  block_id = '1 2'
  block_name = 'domain inclusion'
  second_order = true
[]

[Variables]
  [./disp_x]
    order = SECOND
    family = LAGRANGE
  [../]
  [./disp_y]
    order = SECOND
    family = LAGRANGE
  [../]
  [./pressure]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./mech_x]
    type = LynxADSolidMomentum
    variable = disp_x
    component = 0
  [../]
  [./mech_y]
    type = LynxADSolidMomentum
    variable = disp_y
    component = 1
  [../]
  [./incompressibility]
    type = LynxADMass
    variable = pressure
  [../]
[]

[AuxVariables]
  [./vel_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./vel_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Se]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eqv_strain_rate]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vol_strain_rate]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./vel_x_aux]
    type = LynxVelocityAux
    variable = vel_x
    displacement = disp_x
  [../]
  [./vel_y_aux]
    type = LynxVelocityAux
    variable = vel_y
    displacement = disp_y
  [../]
  [./Se_aux]
    type = LynxADVonMisesStressAux
    variable = Se
  [../]
  [./eqv_strain_rate_aux]
    type = LynxADEqvStrainRateAux
    variable = eqv_strain_rate
  [../]
  [./vol_strain_rate_aux]
    type = LynxADVolStrainRateAux
    variable = vol_strain_rate
  [../]
[]

[Functions]
  [./vel_bc_x]
    type = ParsedFunction
    value = 'vel*x'
    vars = 'vel'
    vals = -0.1
  [../]
  [./vel_bc_y]
    type = ParsedFunction
    value = 'vel*y'
    vars = 'vel'
    vals = 0.1
  [../]
[]

[BCs]
  [./vs_x]
    type = LynxVelocityBC
    variable = disp_x
    boundary = 'left right'
    function = 'vel_bc_x'
  [../]
  [./vs_y]
    type = LynxVelocityBC
    variable = disp_y
    boundary = 'bottom top'
    function = 'vel_bc_y'
  [../]
[]

[Materials]
  [./viscous_domain]
    type = LynxADStokeDeformation
    displacements = 'disp_x disp_y'
    dynamic_pressure = pressure
    block = domain
    A_diffusion = 5.0e-01
  [../]
  [./viscous_inclusion]
    type = LynxADStokeDeformation
    displacements = 'disp_x disp_y'
    dynamic_pressure = pressure
    block = inclusion
    A_diffusion = 5.0e-05
  [../]
[]

[Preconditioning]
  active = 'precond'
  [./precond]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -ksp_max_it -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres asm 1e-10 1e-09 20 50 lu NONZERO'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  start_time = 0
  end_time = 1.0
  dt = 1.0
[]

[Outputs]
  print_linear_residuals = true
  perf_graph = true
  exodus = true
[]
