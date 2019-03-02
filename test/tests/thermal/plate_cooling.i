[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 100
  xmin = 0
  xmax = 1000
  ymin = 0
  ymax = 1000
  zmin = -100000
  zmax = 0
[]

[Variables]
  [./temp]
    initial_condition = 1300
  [../]
[]

[Kernels]
  [./T_time]
    type = TimeDerivative
    variable = temp
  [../]
  [./T_diff]
    type = LynxHeatConduction
    variable = temp
  [../]
[]

[BCs]
  [./T_surface]
    type = DirichletBC
    variable = temp
    boundary = front
    value = 0
  [../]
  [./T_mantle]
    type = DirichletBC
    variable = temp
    boundary = back
    value = 1300
  [../]
[]

[Materials]
  [./density]
    type = LynxDensityConstant
    solid_density = 1000
  [../]
  [./thermo]
    type = LynxThermalConstant
    temperature = temp
    solid_thermal_conductivity = 1
    solid_heat_capacity = 1000
  [../]
[]

[VectorPostprocessors]
  [./line_T]
    type = LineValueSampler
    variable = temp
    start_point = '500 500 0'
    end_point = '500 500 -100000'
    num_points = 101
    sort_by = z
    outputs = csv
  [../]
[]

[Preconditioning]
  active = 'precond'
  [./precond]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -ksp_max_it -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres asm 1E-08 1E-08 20 50 lu NONZERO'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0.0
  end_time = 4.7304e+14 # 15 Myrs
  dt = 1.5768e+13 # 0.5 Myrs
[]

[Outputs]
  print_linear_residuals = true
  perf_graph = true
  exodus = true
  [./csv]
    type = CSV
    sync_times = '3.1536e+13 1.5768e+14 4.7304e+14'
    sync_only = true
  [../]
[]