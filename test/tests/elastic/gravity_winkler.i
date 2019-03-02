[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 2
  ny = 2
  nz = 10
  xmin = -3
  xmax = 3
  ymin = -3
  ymax = 3
  zmin = -30
  zmax = 0
[]

[Functions]
  [./phase_1_fct]
    type = ParsedFunction
    value = 'if(z>=z1, 1.0, 0.0)'
    vars = 'z1'
    vals = '-6'
  [../]
  [./phase_2_fct]
    type = ParsedFunction
    value = 'if(z<z1, if(z>=z2, 1.0, 0.0), 0.0)'
    vars = 'z1 z2'
    vals = '-6 -15'
  [../]
  [./phase_3_fct]
    type = ParsedFunction
    value = 'if(z<z2, 1.0, 0.0)'
    vars = 'z2'
    vals = '-15'
  [../]
  [./pressure_calc_fct]
    type = ParsedFunction
    value = 'if(z>=z1,-rho1*g*z,if(z>=z2,-rho1*g*z1-rho2*g*(z-z1), -rho1*g*z1-rho2*g*(z2-z1)-rho3*g*(z-z2)))'
    vars = 'z1 z2 rho1 rho2 rho3 g'
    vals = '-6 -15 2.5e+03 5.0e+03 8.0e+03 10.0'
  [../]
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
  [./grav_z]
    type = LynxGravityKernel
    variable = disp_z
    component = 2
  [../]
[]

[AuxVariables]
  [./phase_1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      variable = phase_1
      function = phase_1_fct
    [../]
  [../]
  [./phase_2]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      variable = phase_2
      function = phase_2_fct
    [../]
  [../]
  [./phase_3]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      variable = phase_3
      function = phase_3_fct
    [../]
  [../]
  [./pressure_calc]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pressure]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./mises_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_vol]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_dev]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./density]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./pressure_calc_aux]
    type = FunctionAux
    variable = pressure_calc
    function = pressure_calc_fct
    execute_on = 'TIMESTEP_END'
  [../]
  [./pressure_aux]
    type = LynxEffectivePressureAux
    variable = pressure
    execute_on = 'TIMESTEP_END'
  [../]
  [./mises_stress_aux]
    type = LynxVonMisesStressAux
    variable = mises_stress
    execute_on = 'TIMESTEP_END'
  [../]
  [./stress_xx_aux]
    type = LynxStressAux
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = 'TIMESTEP_END'
  [../]
  [./stress_yy_aux]
    type = LynxStressAux
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = 'TIMESTEP_END'
  [../]
  [./stress_zz_aux]
    type = LynxStressAux
    variable = stress_zz
    index_i = 0
    index_j = 0
    execute_on = 'TIMESTEP_END'
  [../]
  [./strain_vol_aux]
    type = LynxVolStrainAux
    variable = strain_vol
    strain_type = 'total'
    execute_on = 'TIMESTEP_END'
  [../]
  [./strain_dev_aux]
    type = LynxEqvStrainAux
    variable = strain_dev
    strain_type = 'total'
    execute_on = 'TIMESTEP_END'
  [../]
  [./density_aux]
    type = MaterialRealAux
    variable = density
    property = bulk_density
    execute_on = 'TIMESTEP_END'
  [../]
[]

[BCs]
  # [./no_x]
  #   type = PresetBC
  #   variable = disp_x
  #   boundary = 'left right'
  #   value = 0.0
  # [../]
  # [./no_y]
  #   type = PresetBC
  #   variable = disp_y
  #   boundary = 'bottom top'
  #   value = 0.0
  # [../]
  [./LynxPressure]
    [./lith_pressure]
      displacements = 'disp_x disp_y disp_z'
      boundary = 'left right bottom top back front'
      function = pressure_calc_fct
    [../]
  [../]
  [./no_z]
    type = PresetBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]
[]

[Materials]
  [./elastic_mat]
    type = LynxDeformation
    displacements = 'disp_x disp_y disp_z'
    compositional_phases = 'phase_1 phase_2 phase_3'
    bulk_modulus = '6.6666667e+09 6.6666667e+09 6.6666667e+09'
    shear_modulus = '4.0e+09 4.0e+09 4.0e+09'
  [../]
  [./density]
    type = LynxDensityConstant
    compositional_phases = 'phase_1 phase_2 phase_3'
    solid_density = '2.5e+03 5.0e+03 8.0e+03'
    gravity_acceleration = 10.0
    has_gravity = true
  [../]
[]

[VectorPostprocessors]
  [./line]
    type = LineValueSampler
    variable = 'pressure pressure_calc'
    num_points = 10
    start_point = '0 0 -1.5'
    end_point = '0 0 -28.5'
    sort_by = 'z'
    outputs = csv
  [../]
[]

[Preconditioning]
  [./hypre]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-snes_linesearch_type -snes_linesearch_maxstep -sneslinesearch_minlambda
                           -snes_atol -snes_rtol -snes_max_it
                           -pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold
                           -pc_hypre_boomeramg_agg_nl
                           -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels
                           -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type
                           -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfactor
                           -ksp_type -ksp_atol -ksp_rtol -ksp_stol -ksp_max_it -ksp_gmres_restart'
    petsc_options_value = 'bt 2e12 1e-3
                           1e-10 1e-20 200
                           hypre boomeramg 0.7
                           4
                           5 25
                           HMIS ext+i
                           2 0.3
                           fgmres 1e-12 1e-15 1e-8 500 30'
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
  execute_on = 'TIMESTEP_END'
  print_linear_residuals = true
  perf_graph = true
  exodus = true
  csv = true
[]
