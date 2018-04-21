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
#ifndef TIME_STEP_CALCULATOR_HPP
#define TIME_STEP_CALCULATOR_HPP

#include "../include/Config.hpp"
#include "../include/Gas.hpp"
#include "../include/Mesh.hpp"
#include "../include/SolutionSnapshot.hpp"

class TimeStepCalculator {
  
  protected:
  
    Config *config_;
    Gas *gas_;
    Mesh *mesh_;
    
    void initializeBase_(Config *config, Gas *gas, Mesh *mesh);
    
  public:
    
    TimeStepCalculator();
    virtual ~TimeStepCalculator();
    
    virtual void initialize(Config *config, Gas *gas, Mesh *mesh) = 0;
    virtual void computeAllTimeStep(
      SolutionSnapshot *solutionSnapshot, const double &courantNum) = 0;
    virtual double getTimeStep(const unsigned long &iElement) = 0;
    
};

#endif // TIME_STEP_CALCULATOR_HPP