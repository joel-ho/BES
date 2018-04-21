################################################################################
# nacaAirfoilGeoGenerator.py
#   Creator: Ho Mun Onn, Joel
#
# Creates .geo file of NACA 4 digit airfoil in rectangular domain.
# Written for python3
################################################################################

import numpy as np

org_coeff = np.asarray((0.2969, -0.126, -0.3516, 0.2843, -0.1015))
mod_coeff = np.asarray((0.2969, -0.126, -0.3516, 0.2843, -0.1036))

def naca_4digit(tc_ratio, coeff, n):
  powers = np.asarray((0.5, 1, 2, 3, 4))
  x = np.linspace(0, 1, n)**2
  y = np.zeros((len(powers), len(x)))
  for i, (c, p) in enumerate(zip(coeff, powers)):
    y[i, :] = c*x**p;
  return x, 5*tc_ratio*np.sum(y, axis=0)

def rot_cw(theta_deg, x, y):
  theta = theta_deg*np.pi/180
  x_new = np.cos(theta)*x + np.sin(theta)*y
  y_new = -np.sin(theta)*x + np.cos(theta)*y
  return x_new, y_new
  
if __name__ == "__main__":

  tunnel_top_left = [-20, 20];
  tunnel_bot_right = [21, -20];
  
  x_air, y_air = naca_4digit(0.12, mod_coeff, 100)
  x_air = np.concatenate((x_air[::-1], x_air[1:-1]))
  y_air = np.concatenate((y_air[::-1], -y_air[1:-1]))
  
  x_air, y_air = rot_cw(0, x_air, y_air)
  
  air_mesh_size = 1/40
  control_volume_edge_mesh_size = 0.75
  overall_mesh_smoothing_radius = 2.5
  boundary_mesh_size = 1.5
  leading_trailing_edge_size_factor = 2.5
  
  with open("nacaAirfoil.geo", "w") as f:
    f.write("////////////////////////////////////////////////////////////////////////////////\n")
    f.write("// nacaAirfoil.geo\n//\tCreator: Ho Mun Onn, Joel\n//\n// NACA 4 digit airfoil in rectangular domain.\n")
    f.write("////////////////////////////////////////////////////////////////////////////////\n\n")
    f.write("Point (4) = {%e, %e, 0, %e};\n"% (tunnel_top_left[0], tunnel_top_left[1], boundary_mesh_size))
    f.write("Point (3) = {%e, %e, 0, %e};\n"% (tunnel_bot_right[0], tunnel_top_left[1], boundary_mesh_size))
    f.write("Point (2) = {%e, %e, 0, %e};\n"% (tunnel_bot_right[0], tunnel_bot_right[1], boundary_mesh_size))
    f.write("Point (1) = {%e, %e, 0, %e};\n"% (tunnel_top_left[0], tunnel_bot_right[1], boundary_mesh_size))
    
    for i in range(1, 4):
      f.write("Line (%d) = {%d, %d};\n"% (i, i, i+1))
    f.write("Line (4) = {4, 1};\n")
  
    i_air_start = 5
    i_air_end = i_air_start
    for i, (x, y) in enumerate(zip(x_air, y_air)):
      f.write("Point (%d) = {%e, %e, 0, %e};\n"% (i_air_start+i, x, y, air_mesh_size))
      i_air_end += 1
    
    air_line_str = "Spline (5) = {"
    for i in range(i_air_start, i_air_end):
      air_line_str = "".join((air_line_str, str(i), ", "))
    air_line_str = "".join((air_line_str, str(i_air_start), "};\n"))
    f.write(air_line_str)
    
    f.write("Line Loop (1) = {1, 2, 3, 4};\n")
    f.write("Line Loop (2) = {5};\n")
    f.write("Plane Surface (1) = {1, 2};\n")

    f.write('Field[1] = Attractor;\n')
    f.write('Field[1].EdgesList = {5};\n')

    f.write('Field[2] = MathEval;\n')
    f.write('Field[2].F = \"%f*(F1/%f)^2 + %f\";\n'% (control_volume_edge_mesh_size-air_mesh_size, overall_mesh_smoothing_radius, air_mesh_size))

    f.write('Field[3] = Attractor;\n')
    f.write('Field[3].NodesList = {%d, %d};\n'% (i_air_start, np.round((i_air_end-i_air_start)/2)))
    
    f.write('Field[4] = MathEval;\n')
    f.write('Field[4].F = \"%f*(F3/%f)^2 + %f\";\n'% (air_mesh_size-air_mesh_size/leading_trailing_edge_size_factor, 0.2, air_mesh_size/leading_trailing_edge_size_factor))
    
    f.write('Field[5] = Min;\n')
    f.write('Field[5].FieldsList = {2, 4};\n')
    
    f.write('Background Field = 5;\n')
    
    f.write("Physical Line (\"FREESTREAM\") = {1, 2, 3, 4};\n")
    f.write("Physical Line (\"SLIP_WALL\") = {5};\n")
    
    f.write("Physical Surface (\"FLUID\") = {1};\n")
    