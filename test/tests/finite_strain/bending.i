[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 10
  ymin = 0
  ymax = 2
  nx = 10
  ny = 2
  elem_type = QUAD4
[]

[MeshModifiers]
  [./corner]
    type = AddExtraNodeset
    new_boundary = 101
    coord = '0 0'
  [../]
  [./side]
    type = AddExtraNodeset
    new_boundary = 102
    coord = '10 0'
  [../]
  [./mid]
    type = AddExtraNodeset
    new_boundary = 103
    coord = '5 2'
  [../]
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
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

[BCs]
 [./fix_corner_x]
   type = PresetBC
   variable = disp_x
   boundary = 101
   value = 0
 [../]
 [./fix_corner_y]
   type = PresetBC
   variable = disp_y
   boundary = 101
   value = 0
 [../]
 [./fix_y]
   type = PresetBC
   variable = disp_y
   boundary = 102
   value = 0
 [../]
 [./move_y]
   type = FunctionPresetBC
   variable = disp_y
   boundary = 103
   function = '-t'
 [../]
[]

[Materials]
  [./elastic_mat]
    type = LynxADElasticDeformation
    displacements = 'disp_x disp_y'
    bulk_modulus = 6.7866666667e+04
    shear_modulus = 0.754e+05
    strain_model = finite
  [../]
[]

[Preconditioning]
  [./precond]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type -snes_atol'
    petsc_options_value = 'hypre boomeramg 1.0e-12'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  dt = 0.1
  dtmin = 0.1
  num_steps = 2
[]

[Outputs]
  print_linear_residuals = true
  perf_graph = true
  exodus = true
[]
