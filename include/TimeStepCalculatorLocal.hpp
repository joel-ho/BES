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
#ifndef TIME_STEP_CALCULATOR_LOCAL_HPP
#define TIME_STEP_CALCULATOR_LOCAL_HPP

#include <cmath>

#include "../include/TimeStepCalculator.hpp"
#include "../include/GasPerfect.hpp"

class TimeStepCalculatorLocal: public TimeStepCalculator {
  
  protected:
  
    double *allTimeStep_;
    double **projectedLength_;
    
    void computeProjectedLength_();
    
  public:
    ~TimeStepCalculatorLocal();
    void initialize(Config *config, Gas *gas, Mesh *mesh);
    void computeAllTimeStep(
      SolutionSnapshot *solutionSnapshot, const double &courantNum);
    double getTimeStep(const unsigned long &iElement);
    
};

#endif // TIME_STEP_CALCULATOR_LOCAL_HPP