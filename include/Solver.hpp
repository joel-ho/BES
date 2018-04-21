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
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <map>

#include "../include/Options.hpp"
#include "../include/Config.hpp"
#include "../include/Gas.hpp"
#include "../include/Mesh.hpp"
#include "../include/Solution.hpp"
#include "../include/Gradient.hpp"
#include "../include/Limiter.hpp"
#include "../include/Flux.hpp"
#include "../include/FluxJacobian.hpp"
#include "../include/TimeStepCalculator.hpp"

class Solver {
  
  protected:
    
    Config *config_;
    Gas *gas_;
    Mesh *mesh_;
    Solution *solution_;
    
    Gradient *gradient_;
    Limiter *limiter_;
    TimeStepCalculator *timeStepCalculator_;
    
    void initializeBase_(Config *config, Gas *gas, Mesh *mesh);
    void getFirstOrderFaceVariable_(
      SolutionSnapshot *solutionSnapshot,
      unsigned long &iFace,
      double * leftVariable, 
      double * rightVariable);
    void computeFaceConservedVariable_(
      SolutionSnapshot *solutionSnapshot,
      Gradient *gradient,
      Limiter *limiter,
      unsigned long &iFace,
      double *leftVariable, 
      double *rightVariable);
    void writeToFile_(const unsigned long &iTimeStep, const short &iCurrentSnapshot);
    void appendToResidual_(string residualFileName, string residual);
    
  public:
  
    Solver();
    virtual ~Solver();
    
    virtual void initialize(Config *config, Gas *gas, Mesh *mesh) = 0;
    virtual void setUniversalInitialCondition(double *primitiveVariable) = 0;
    virtual void solve() = 0;
    
};

#endif // SOLVER_HPP