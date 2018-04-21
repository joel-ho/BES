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
#ifndef SOLUTION_SNAPSHOT_HPP
#define SOLUTION_SNAPSHOT_HPP

#include <fstream>
#include <iostream>

#include "../include/Options.hpp"
#include "../include/Gas.hpp"
#include "../include/GasPerfect.hpp"
#include "../include/Config.hpp"
#include "../include/Mesh.hpp"

using namespace std;

class SolutionSnapshot {
  
  private:
    
    double ***variableContainer_;
    double **conservedVariables_, **primitiveVariables_;
    double **dimensionalPrimitiveVariables_;
    
    Config *config_;
    Gas *gas_;
    Mesh *mesh_;
    
  public:
    
    SolutionSnapshot();
    ~SolutionSnapshot();
    
    void initialize(Config *config, Mesh *mesh, Gas *gas);
    
    void setConservedVariable(
      const unsigned long &iElement, const short &iVariable, double value);
    void setConservedVariable (
      const unsigned long &iElement, double *conservedVariable);
    
    void computeAllPrimitiveVariable();
    
    double* getConservedVariable(const unsigned long &iElement);
    double getConservedVariable(
      const unsigned long &iElement, const short &iVariable);
    double getConservedVariable(
      const unsigned long &iElement, FlowVariable flowVariable);
      
    void addToConservedVariable(
      const unsigned long &iElement, const short &iVariable, double value);
      
    // No check on whether primitive variables are initialized
    double* getPrimitiveVariable(const unsigned long &iElement);
    double getPrimitiveVariable(
      const unsigned long &iElement, const short &iVariable);
    double getPrimitiveVariable(
      const unsigned long &iElement, FlowVariable flowVariable);
      
    double getVariable(
      const unsigned long &iElement, const short &iVariable, VariableType variableType);
      
    void writePrimitiveToCsv(string fileName);
    void writePrimitiveToVtk(string fileName);
    void writeSlipWallToCsv(string fileName);
    
};

#endif // SOLUTION_SNAPSHOT_HPP