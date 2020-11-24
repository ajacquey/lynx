[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 40
  nx = 160
[]

[Functions]
  [./initial_phase]
    type = ParsedFunction
    value = 'if(x>=5&x<=10,1,if(x>=15&x<=20,if(x<=17.5,(x-15)/2.5,1.0-(x-17.5)/2.5),0))'
  [../]
[]

[Variables]
  [./phase]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      variable = phase
      function = initial_phase
    [../]
  [../]
[]

[LynxAdvection]
  velocities = 'vel_x'
  compositional_phases = 'phase'
  beta_stabilization = 0.026
  cr_stabilization = 1
[]

[BCs]
  [./left_phase]
    type = DirichletBC
    variable = phase
    value = 0.0
    boundary = 'left right'
    preset = true
  [../]
[]

[AuxVariables]
  [./vel_x]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1.0
  [../]
  [./vel_norm]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  [./vel_norm_aux]
    type = LynxVelocityNormAux
    variable = vel_norm
    velocities = vel_x
    execute_on = 'timestep_end'
  [../]
[]

[Postprocessors]
  [./cfl]
    type = LynxExplicitTimeStepSelector
    beta = 0.1
    vel_norm = vel_norm
    has_premult = false
    initial_value = 1e-5
    maximum_value = 1
    execute_on = 'initial timestep_end'
  [../]
[]

[Preconditioning]
  active = 'tight_asm'
  [./loose_asm]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-snes_linesearch_type -snes_atol -snes_rtol -snes_max_it
                            -ksp_type -ksp_atol -ksp_max_it
                            -pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'basic 1e-12 1e-15 200
                           fgmres 1e-4 500
                           asm ilu NONZERO'
  [../]
  [./tight_asm]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname =  '-snes_atol -snes_rtol -snes_max_it
                           -ksp_type -ksp_atol -ksp_max_it
                           -pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = '1e-10 1e-20 200
                           fgmres 1e-8 500
                           asm ilu NONZERO'
  [../]
  [./tight_hypre]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-snes_atol -snes_rtol -snes_max_it
                           -ksp_type -ksp_atol -ksp_max_it
                           -pc_type -pc_hypre_type
                           -pc_hypre_boomeramg_coarsen_type
                           -pc_hypre_boomeramg_agg_nl'
    petsc_options_value = '1e-16 1e-20 200
                           fgmres 1e-8 150
                           hypre boomeramg
                           HMIS
                           2'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = bdf2
  start_time = 0
  end_time = 15
  [./TimeStepper]
    type = PostprocessorDT
    postprocessor = cfl
  [../]
[]

[Outputs]
  print_linear_residuals = false
  perf_graph = true
  [./console]
    type = Console
    execute_postprocessors_on = none
  [../]
  [./out]
    type = Exodus
    interval = 10
  [../]
[]
