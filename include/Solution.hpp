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
#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <sstream>
#include <map>

#include "../include/Config.hpp"
#include "../include/Gas.hpp"
#include "../include/GasPerfect.hpp"
#include "../include/Mesh.hpp"
#include "../include/SolutionSnapshot.hpp"

using namespace std;

class Solution {
  
  private:
    
    short nCurrentSnapshots_;
    
    SolutionSnapshot **currentSnapshots_;
    
    Config* config_;
    Gas *gas_;
    Mesh* mesh_;
    
    short indexCurrentSnapshot_(
      const unsigned long &iTimeStep, const short &iCurrentSnapshot);
    
  public:
    
    Solution();
    ~Solution();
    
    void initialize(Config* config, Gas *gas, Mesh* mesh);
    
    void setConservedVariable(
      const unsigned long &iTimeStep, const short &iCurrentSnapshot,
      const unsigned long &iElement, double *conservedVariable);
    void setConservedVariable(
      const unsigned long &iTimeStep, const short &iCurrentSnapshot,
      const unsigned long &iElement, const short &iVariable,
      double value);
    double* getConservedVariable(
      const unsigned long &iTimeStep, const short &iCurrentSnapshot, 
      const unsigned long &iElement);
    double getConservedVariable(
      const unsigned long &iTimeStep, const short &iCurrentSnapshot, 
      const unsigned long &iElement, const short &iVariable);
      
    // No check on whether primitive variables are initialized
    double getPrimitiveVariable(
      const unsigned long &iTimeStep, const short &iCurrentSnapshot, 
      const unsigned long &iElement, const short &iVariable);
      
    
    SolutionSnapshot* getCurrentSnapshot(
      const unsigned long &iTimeStep, const short &iCurrentSnapshot);
    void writeSolution(const unsigned long &iTimeStep, const short &iCurrentSnapshot);
    void writeSlipWall(const unsigned long &iTimeStep, const short &iCurrentSnapshot);
    
};

#endif // SOLUTION_HPP