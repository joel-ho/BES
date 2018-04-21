################################################################################
# forwardStepGeoGenerator.py
#   Creator: Ho Mun Onn, Joel
#
# Creates .geo file of Woodward and Colella forward step problem with radius
# at singularity point.
# Written for python3.
################################################################################

import numpy as np

mesh_size = 1/60
tip_mesh_size = mesh_size/8

tip_fillet = 1/80
tip_refinement_zone_radius = 0.2
n_fillet_pts = 30

with open("forwardStep.geo", "w") as f:
  
  f.write("////////////////////////////////////////////////////////////////////////////////\n")
  f.write("// forwardStep.geo\n//\tCreator: Ho Mun Onn, Joel\n//\n// Woodward and Colella forward step problem.\n")
  f.write("////////////////////////////////////////////////////////////////////////////////\n\n")
  f.write("Point (1) = {0, 0, 0, %e};\n"% (mesh_size))
  f.write("Point (2) = {0.6, 0, 0, %e};\n"% (mesh_size))
  f.write("Point (3) = {0.6, %e, 0, %e};\n"% (0.2-tip_fillet, mesh_size))
  f.write("Point (4) = {%e, 0.2, 0, %e};\n"% (0.6+tip_fillet, mesh_size))
  f.write("Point (5) = {3, 0.2, 0, %e};\n"% (mesh_size))
  f.write("Point (6) = {3, 1, 0, %e};\n"% (mesh_size))
  f.write("Point (7) = {0, 1, 0, %e};\n"% (mesh_size))

  for i in range(1, 7):
    if i != 3:
      f.write("Line (%d) = {%d, %d};\n"% (i, i, i+1))
  f.write("Line (7) = {7, 1};\n")
  
  angle = np.linspace(np.pi/2, np.pi, n_fillet_pts+2)
  angle = angle[1:-1][::-1]
  for i in range(0, n_fillet_pts):
    circle_x = tip_fillet*np.cos(angle[i]) + 0.6+tip_fillet
    circle_y = tip_fillet*np.sin(angle[i]) + 0.2-tip_fillet
    f.write("Point (%d) = {%e, %e, 0, %e};\n"% (8+i, circle_x, circle_y, mesh_size))
  
  f.write("Spline (8) = {3, ")
  for i in range(0, n_fillet_pts):
    f.write("%d, "% (8+i))
  f.write("4};")
  
  f.write("Line Loop (1) = {1, 2, 8, 4, 5, 6, 7};\n")
  f.write("Plane Surface (1) = {1};\n")
  
  f.write('Field[1] = Attractor;\n')
  f.write('Field[1].NodesList = {3};\n')

  f.write('Field[2] = MathEval;\n')
  f.write('Field[2].F = \"%f*(F1/%f)^2 + %f\";\n'% (mesh_size-tip_mesh_size, tip_refinement_zone_radius, tip_mesh_size))
  
  f.write('Background Field = 2;\n')
  
  f.write("Physical Line (\"SLIP_WALL\") = {1, 2, 8, 4, 6};\n")
  f.write("Physical Line (\"FREESTREAM\") = {5, 7};\n")
  
  f.write("Physical Surface (\"FLUID\") = {1};\n")