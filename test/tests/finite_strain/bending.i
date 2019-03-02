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

[GlobalParams]
  displacements = 'disp_x disp_y'
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
    type = LynxSolidMomentum
    variable = disp_x
    component = 0
  [../]
  [./mech_y]
    type = LynxSolidMomentum
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
    type = LynxDeformation
    bulk_modulus = 6.7866666667e+04
    shear_modulus = 0.754e+05
    strain_model = finite
  [../]
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '_pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-10
  nl_max_its = 10

  l_tol  = 1e-4
  l_max_its = 50

  dt = 0.1
  dtmin = 0.1

  num_steps = 2
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Outputs]
  exodus = true
[]
