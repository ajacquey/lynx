[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 4
  nz = 1
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  zmin = 0
  zmax = 0.25
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
  [./damage]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mises_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eqv_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eqv_in_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eqv_in_strain_rate]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./damage_yield_function]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./damage_aux]
    type = ConstantAux
    variable = damage
    value = 0.0
  [../]
  [./mises_stress_aux]
    type = LynxVonMisesStressAux
    variable = mises_stress
    execute_on = 'timestep_end'
  [../]
  [./eqv_strain_aux]
    type = LynxEqvStrainAux
    variable = eqv_strain
    execute_on = 'timestep_end'
  [../]
  [./eqv_in_strain_aux]
    type = LynxEqvStrainAux
    variable = eqv_in_strain
    strain_type = inelastic
    execute_on = 'timestep_end'
  [../]
  [./eqv_in_strain_rate_aux]
    type = LynxEqvStrainRateAux
    variable = eqv_in_strain_rate
    strain_type = inelastic
    execute_on = 'timestep_end'
  [../]
  [./damage_yield_function_aux]
    type = MaterialRealAux
    variable = damage_yield_function
    property = damage_yield_function
  [../]
[]

[BCs]
  [./no_ux]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./ux_right]
    type = LynxVelocityBC
    variable = disp_x
    boundary = right
    value = -1.0e-14
  [../]
  [./uy_top]
    type = LynxVelocityBC
    variable = disp_y
    boundary = top
    value = 1.0e-14
  [../]
  [./no_uy]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./no_uz]
    type = PresetBC
    variable = disp_z
    boundary = 'front back'
    value = 0.0
  [../]
[]

[Materials]
  [./damage_deformation]
    type = LynxDamageDeformation
    displacements = 'disp_x disp_y disp_z'
    damage = damage
    bulk_modulus = 1.0e+10
    shear_modulus = 1.0e+10
    damage_modulus = 0.0
    friction_angle = 30
    cohesion = 10.0e+06
    plastic_viscosity = 1.0e+21
  [../]
[]

[Postprocessors]
  [./Se]
    type = ElementAverageValue
    variable = mises_stress
    outputs = csv
  [../]
  [./E_eqv]
    type = ElementAverageValue
    variable = eqv_strain
    outputs = csv
  [../]
  [./E_in_eqv]
    type = ElementAverageValue
    variable = eqv_in_strain
    outputs = csv
  [../]
  [./E_dot_in_eqv]
    type = ElementAverageValue
    variable = eqv_in_strain_rate
    outputs = csv
  [../]
  [./Yield]
    type = ElementAverageValue
    variable = damage_yield_function
    outputs = csv
  [../]
[]

[Preconditioning]
  [./precond]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor -snes_atol'
    petsc_options_value = 'hypre boomeramg 0.7 4 5 25 HMIS ext+i 2 0.3 1.0e-08'
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