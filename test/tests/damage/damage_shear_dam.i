[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 5
  ny = 5
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
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
[]

[AuxVariables]
  [./p_lith]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1.0e+01
  [../]
  [./mises_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vol_strain]
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
  [./yield]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./damage]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./damage_force]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pressure]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./p_lith_aux]
    type = ConstantAux
    variable = p_lith
    value = 1.0e+01
  [../]
  [./mises_stress_aux]
    type = LynxADVonMisesStressAux
    variable = mises_stress
    execute_on = 'TIMESTEP_END'
  [../]
  [./vol_strain_aux]
    type = LynxADVolStrainAux
    variable = vol_strain
    execute_on = 'TIMESTEP_END'
  [../]
  [./eqv_strain_aux]
    type = LynxADEqvStrainAux
    variable = eqv_strain
    execute_on = 'TIMESTEP_END'
  [../]
  [./eqv_in_strain_aux]
    type = LynxADEqvStrainAux
    variable = eqv_in_strain
    strain_type = plastic
    execute_on = 'TIMESTEP_END'
  [../]
  [./eqv_in_strain_rate_aux]
    type = LynxADEqvStrainRateAux
    variable = eqv_in_strain_rate
    strain_type = plastic
    execute_on = 'TIMESTEP_END'
  [../]
  [./yield_aux]
    type = ADMaterialRealAux
    variable = yield
    property = plastic_yield_function
    execute_on = 'TIMESTEP_END'
  [../]
  [./damage_aux]
    type = ADMaterialRealAux
    variable = damage
    property = damage
    execute_on = 'TIMESTEP_END'
  [../]
  [./damage_force_aux]
    type = ADMaterialRealAux
    variable = damage_force
    property = damage_force
    execute_on = 'TIMESTEP_END'
  [../]
  [./pressure_aux]
    type = LynxADEffectivePressureAux
    variable = pressure
    execute_on = 'TIMESTEP_END'
  [../]
[]

[BCs]
  [./no_ux]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
    preset = true
  [../]
  [./ux_right]
    type = LynxVelocityBC
    variable = disp_x
    boundary = right
    value = -1.0
  [../]
  [./uy_top]
    type = LynxVelocityBC
    variable = disp_y
    boundary = top
    value = 1.0
  [../]
  [./no_uy]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
    preset = true
  [../]
[]


[Materials]
  [./elastic_mat]
    type = LynxADDamageDeformation
    displacements = 'disp_x disp_y'
    lithostatic_pressure = 'p_lith'
    bulk_modulus = 1.0e+04
    shear_modulus = 1.0e+04
    damage_model = 'damage_vladi'
  [../]
  [./damage_vladi]
    type = LynxADLyakhovskyDamage
    damage_modulus = 1.0e+04
    friction_angle = 30
    plastic_viscosity = 1.0e-01
    initial_damage = 1.0e-03
  [../]
[]

[Postprocessors]
   [./D]
    type = ElementAverageValue
    variable = damage
    outputs = csv
  [../]
  [./Se]
    type = ElementAverageValue
    variable = mises_stress
    outputs = csv
  [../]
  [./E_vol]
    type = ElementAverageValue
    variable = vol_strain
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
  [./F]
    type = ElementAverageValue
    variable = yield
    outputs = csv
  [../]
  [./Dam_force]
    type = ElementAverageValue
    variable = damage_force
    outputs = csv
  [../]
[]

[Preconditioning]
  [./precond]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_type -snes_atol -snes_rtol -snes_max_it -ksp_max_it
                           -pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold
                           -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths
                           -pc_hypre_boomeramg_max_levels
                           -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type
                           -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor
                           -snes_linesearch_type'
    petsc_options_value = 'fgmres 1.0e-06 1.0e-10 1000 100
                           hypre boomeramg 0.7
                           4 5
                           25
                           HMIS ext+i
                           2 0.3
                           basic'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0.0
  end_time = 3.0e-03
  num_steps = 50
[]

[Outputs]
  execute_on = 'timestep_end'
  print_linear_residuals = false
  perf_graph = true
  exodus = true
  csv = true
[]
