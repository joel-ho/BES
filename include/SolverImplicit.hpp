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
#ifndef SOLVER_IMPLICIT_HPP
#define SOLVER_IMPLICIT_HPP

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
#include "../include/FluxJacobianWrapper.hpp"
#include "../include/TimeStepCalculatorLocal.hpp"
#include "../include/TimeStepCalculatorGlobal.hpp"
#include "../include/CoefficientMatrix.hpp"

struct FaceLoopVariables {
  // Private variables in parfor loop. New instance created via copy 
  // constructor through firstprivate directive when 
  // parfor loop is called.
  
  // Individual buffers and flux calculators needed for each thread
  unsigned long leftElementId, rightElementId, iTimeStep;
  double faceArea;
  double leftVariableBuffer[VAR_BUFFER];
  double rightVariableBuffer[VAR_BUFFER];
  double fluxBuffer[VAR_BUFFER];
  double leftJacobianBuffer[VAR_BUFFER*VAR_BUFFER];
  double rightJacobianBuffer[VAR_BUFFER*VAR_BUFFER];
  
  map <BoundaryType, Flux*> fluxCalculator;
  map <BoundaryType, FluxJacobian*> fluxJacobianCalculator;
  map <BoundaryType, double*> rightVariableMap;
  map <BoundaryType, double*> rightJacobianMap;
  
  BoundaryType boundaryType;
  Face *face;
  Element *leftElement, *rightElement;
  
  // Variables to aid initialization
  Config *config;
  Gas *gas;
  Mesh *mesh;
  
  // Functions
  FaceLoopVariables();
  ~FaceLoopVariables();
  
  void initialize(Config* configIn, Gas* gasIn, Mesh* meshIn, const unsigned long iTimeStepIn, int printOptions);
  
  FaceLoopVariables(const FaceLoopVariables& solverImplicitVariables);
  
  void printJacobian(short nVar, double* jacobian);
  
};

class SolverImplicit: public Solver {
  
  private:
    
    double actualCourantNum_;
    double * maxResidual_, * maxResidualOld_;
    vector<double> explicitFlux_, deltaConservedVariable_;
    
    CoefficientMatrix* coefficientMatrix_;
    
    FaceLoopVariables faceLoopVariables_;
    
    unsigned long computeVectorIdx_(unsigned long iElement, short iVar);
    
    void timeStep_(
      const unsigned long & iTimeStepStart, 
      FaceLoopVariables loopVar, 
      string residualFileName);
      
    void resetSystemOfEquation_();
      
    void computeFaceFluxes_(
      const unsigned long & iTimeStepStart,
      SolutionSnapshot * solutionSnapshot,
      FaceLoopVariables loopVar);
      
    void updateSolution_(SolutionSnapshot * solutionSnapshot);
    
    void adaptCourantNum_();
    
  public:
    
    ~SolverImplicit();
  
    void initialize(Config *config, Gas *gas, Mesh *mesh);
    void setUniversalInitialCondition(double *primitiveVariable);
    void solve();
    
};

#endif // SOLVER_IMPLICIT_HPP
