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
#ifndef SOLVER_EXPLICIT_HPP
#define SOLVER_EXPLICIT_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>

#include "../include/Solver.hpp"
#include "../include/GasPerfect.hpp"
#include "../include/GradientGreenGauss.hpp"
#include "../include/LimiterVenkatakrishnan.hpp"
#include "../include/FluxRoe.hpp"
#include "../include/FluxFreestream.hpp"
#include "../include/FluxSlipWall.hpp"
#include "../include/TimeStepCalculatorLocal.hpp"
#include "../include/TimeStepCalculatorGlobal.hpp"

class PrivateVariableWrapper {
  
  public:
    
    short nVariable;
    unsigned long iTimeStep;
    double *leftVariable, *rightVariable, *flux;
    map <BoundaryType, Flux*> fluxCalculator;
    map <BoundaryType, double*> boundaryVariable;
    
    Config *config;
    Gas *gas;
    Face *face;
    Element *element;
    
    PrivateVariableWrapper();
    PrivateVariableWrapper(const PrivateVariableWrapper &wrapper);
    ~PrivateVariableWrapper();
    
};

class SolverExplicit: public Solver {
  
  protected:
  
    double **subResidual_;
    double *maxResidual_;
  	
    PrivateVariableWrapper privateVariable_;
  	
    void clearDeltaConservedVariable_();
    void timeStep_(
      const unsigned long &iTimeStepStart, PrivateVariableWrapper privateVariable, string residualFileName);
    
  public:
    
    ~SolverExplicit();
  
    void initialize(Config *config, Gas *gas, Mesh *mesh);
    void setUniversalInitialCondition(double *primitiveVariable);
    void solve();
    
};

#endif // SOLVER_EXPLICIT_HPP
