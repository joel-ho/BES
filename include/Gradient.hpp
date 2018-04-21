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
#ifndef GRADIENT_HPP
#define GRADIENT_HPP

#include "../include/Options.hpp"
#include "../include/Config.hpp"
#include "../include/Mesh.hpp"
#include "../include/SolutionSnapshot.hpp"

class Gradient {
  
  protected:
  
    short nVariable_;
    double ***variableGradient_;
    
    VariableType variableType_;
    
    Config* config_;
    Mesh* mesh_;
    
  public:
    Gradient();
    virtual ~Gradient();
    
    void initialize(Config* config, Mesh* mesh, VariableType variableType);
    
    virtual void computeGradient(SolutionSnapshot *solutionSnapshot) = 0;
    
    double getGradient(
      const unsigned long &iElement, const short &iVariable, const short &iDim);
      
};

#endif // GRADIENT_HPP
