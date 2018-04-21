/*
  Basic Euler Solver (BES) - A lightweight Euler Solver for fluid dynamics.
  Copyright (C) 2018, Joel Ho Mun Onn

  BES is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  BES is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with BES.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <string>
#include <map>

#define LINE_BUFFER 256
#define DIM_BUFFER 3
#define VAR_BUFFER 10

using namespace std;

enum class MeshFormat { GMSH=1 };
  
enum class BoundaryType { 
  INTERNAL=1, SLIP_WALL, FREESTREAM };
  
extern map <string, BoundaryType> BoundaryTypeStringMap;

extern map <BoundaryType, int> ExternalBoundaryMap;

enum class EquationSet { EULER_NON_DIMENSIONAL=1 };

enum class NumVariable { 
  EULER_TWO_DIM_CONSERVATIVE=4, EULER_TWO_DIM_PRIMITIVE=5, 
  EULER_THREE_DIM_CONSERVATIVE=5, EULER_THREE_DIM_PRIMITIVE=6 };

enum class SolverType { EXPLICIT=1, IMPLICIT };

enum class AdaptiveCourantNum { YES=1, NO };

enum class LocalTimeStep { YES=1, NO };

enum class ExplicitScheme { MULTI_STEP=1, EULER_FORWARD };

enum class SpatialDiscretizationOrder { FIRST_ORDER=1, SECOND_ORDER };

enum class GradientScheme { GREEN_GAUSS=1 };

enum class LimiterScheme { VENKATAKRISHNAN=1 };

enum class FluxDiscretizationScheme { ROE=1 };

enum class AmgPrecondCoarsener { AGGREGATION=1, SMOOTHED_AGGREGATION, SMOOTHED_AGGREGATION_EMIN, RUGE_STUBEN };

enum class AmgPrecondSmoother { DAMPED_JACOBI=1, ILU0, SPAI0, GAUSS_SEIDEL };

enum class LinearSolver { FGMRES=1, GMRES, LGMRES, BICGSTAB };

enum class FlowVariable { 
  RHO=0, U, V, W, P, T, P_TOT, T_TOT, 
  RHO_U, RHO_V, RHO_W, RHO_E, H, C, V_CONTRAVARIANT, Q_SQ,
  COUNT};
  
enum VariableType { CONSERVED=0, PRIMITIVE, COUNT };

enum class OutputFormat { VTK=1, CSV };

enum class OutputSlipWall { YES=1, NO };

#endif // OPTIONS_HPP