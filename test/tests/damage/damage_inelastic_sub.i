[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

[Variables]
  [./damage]
  [../]
[]

[Kernels]
  [./damage_time]
    type = TimeDerivative
    variable = damage
  [../]
  [./damage_rate]
    type = LynxDamageRate
    variable = damage
    damage_rate = damage_rate
  [../]
[]

[AuxVariables]
  [./damage_rate]
    order = CONSTANT
    family = MONOMIAL
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
    petsc_options_value = 'fgmres 1.0e-25 1.0e-10 1000 100
                           hypre boomeramg 0.7
                           4 5
                           25
                           HMIS ext+i
                           2 0.3
                           bt'
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
  print_linear_residuals = false
  perf_graph = false
  [./console]
    type = Console
    execute_postprocessors_on = none
  [../]
[]
