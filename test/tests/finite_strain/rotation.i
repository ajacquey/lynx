#
# Rotation Test
#
# This test is designed to compute stress based on uniaxial strain
# and then follow that stress as the mesh is rotated 90 degrees.
#
# The mesh is composed of one block with a single element.  The nodal
# displacements in the three directions are prescribed.  Poisson's
# ratio is 1/3, and Young's modulus is 1e6.
#
# This test is mentioned in
# K. Kamojjala, R. Brannon, A. Sadeghirad, and J. Guilkey, "Verification
#   tests in solid mechanics," Engineering with Computers, Vol. 31, 2015.
#   DOI: 10.1007/s00366-013-0342-x
#
[Mesh]
  type = FileMesh
  file = rotation_test.e
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
    type = LynxADSolidMomentum
    variable = disp_x
    component = 0
  [../]
  [./mech_y]
    type = LynxADSolidMomentum
    variable = disp_y
    component = 1
  [../]
  [./mech_z]
    type = LynxADSolidMomentum
    variable = disp_z
    component = 2
  [../]
[]

[AuxVariables]
  [./Sxx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Syy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Szz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Sxy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Syz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Szx]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./Sxx_aux]
    type = LynxStressAux
    variable = Sxx
    index_i = 0
    index_j = 0
  [../]
  [./Syy_aux]
    type = LynxStressAux
    variable = Syy
    index_i = 1
    index_j = 1
  [../]
  [./Szz_aux]
    type = LynxStressAux
    variable = Szz
    index_i = 2
    index_j = 2
  [../]
  [./Sxy_aux]
    type = LynxStressAux
    variable = Sxy
    index_i = 0
    index_j = 1
  [../]
  [./Syz_aux]
    type = LynxStressAux
    variable = Syz
    index_i = 1
    index_j = 2
  [../]
  [./Szx_aux]
    type = LynxStressAux
    variable = Szx
    index_i = 2
    index_j = 0
  [../]
[]

[Functions]
  [./x_200]
    type = ParsedFunction
    vars = 'delta t0'
    vals = '1e-6 1.0'
    value = 'if(t<=1.0, delta*t, (1.0+delta)*cos(pi/2*(t-t0)) - 1.0)'
  [../]
  [./y_200]
    type = ParsedFunction
    vars = 'delta t0'
    vals = '1e-6 1.0'
    value = 'if(t<=1.0, 0.0, (1.0+delta)*sin(pi/2*(t-t0)))'
  [../]
  [./x_300]
    type = ParsedFunction
    vars = 'delta t0'
    vals = '1e-6 1.0'
    value = 'if(t<=1.0, delta*t, (1.0+delta)*cos(pi/2.0*(t-t0)) - sin(pi/2.0*(t-t0)) - 1.0)'
  [../]
  [./y_300]
    type = ParsedFunction
    vars = 'delta t0'
    vals = '1e-6 1.0'
    value = 'if(t<=1.0, 0.0, cos(pi/2.0*(t-t0)) + (1+delta)*sin(pi/2.0*(t-t0)) - 1.0)'
  [../]
  [./x_400]
    type = ParsedFunction
    vars = 'delta t0'
    vals = '1e-6 1.0'
    value = 'if(t<=1.0, 0.0, -sin(pi/2.0*(t-t0)))'
  [../]
  [./y_400]
    type = ParsedFunction
    vars = 'delta t0'
    vals = '1e-6 1.0'
    value = 'if(t<=1.0, 0.0, cos(pi/2.0*(t-t0)) - 1.0)'
  [../]
[]

[BCs]
  [./no_x]
    type = PresetBC
    variable = disp_x
    boundary = 100
    value = 0.0
  [../]
  [./no_y]
    type = PresetBC
    variable = disp_y
    boundary = 100
    value = 0.0
  [../]
  [./x_200]
    type = FunctionPresetBC
    variable = disp_x
    boundary = 200
    function = x_200
  [../]
  [./y_200]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 200
    function = y_200
  [../]
  [./x_300]
    type = FunctionPresetBC
    variable = disp_x
    boundary = 300
    function = x_300
  [../]
  [./y_300]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 300
    function = y_300
  [../]
  [./x_400]
    type = FunctionPresetBC
    variable = disp_x
    boundary = 400
    function = x_400
  [../]
  [./y_400]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 400
    function = y_400
  [../]
  [./no_z]
    type = PresetBC
    variable = disp_z
    boundary = '100 200 300 400'
    value = 0.0
  [../]
[]

[Materials]
  [./elastic_mat]
    type = LynxADElasticDeformation
    displacements = 'disp_x disp_y disp_z'
    strain_model = finite
    bulk_modulus = 1.0e+06
    shear_modulus = 0.375e+06
  [../]
[]

[Preconditioning]
  [./precond]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -ksp_max_it -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres asm 1E-10 1E-15 20 50 lu NONZERO'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  start_time = 0.0
  dt = 0.01
  end_time = 2.0
[]

[Outputs]
  print_linear_residuals = true
  perf_graph = true
  exodus = true
[]
