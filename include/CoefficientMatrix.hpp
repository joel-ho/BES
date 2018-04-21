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
#ifndef COEFFICIENT_MATRIX_HPP
#define COEFFICIENT_MATRIX_HPP

#include <vector>

#include <amgcl/runtime.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>

#include <amgcl/coarsening/aggregation.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/smoothed_aggr_emin.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>

#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>

#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/lgmres.hpp>
#include <amgcl/solver/fgmres.hpp>
#include <amgcl/solver/bicgstab.hpp>

#include "../include/MatrixSparse.hpp"
#include "../include/Options.hpp"
#include "../include/Config.hpp"
#include "../include/Mesh.hpp"

class CoefficientMatrix {
  
  protected:
    
    unsigned long nBlock_;
    short blockSize_;
    
    MatrixSparse* matrix_;
    
    Config* config_;
    Mesh* mesh_;
    
    void printOptions_();
    
    unsigned long computeIdx_(unsigned long iElement, short iVar);
    
  public:
    CoefficientMatrix ();
    ~CoefficientMatrix ();
    
    void initialize(Config* config, Mesh* mesh);
    
    void changeAllVal(double newVal);
    void addToDiagonal(double value);
    void addToVal(
      unsigned long iElementRow, unsigned long iElementCol, 
      short iVarRow, short iVarCol, 
      double addVal);
    void addToBlock(
      unsigned long iElementRow, unsigned long iElementCol, 
      double multiplier, double * addValRowMajor);
    
    MatrixSparse* getMatrixSparse();
    
    void solve(vector<double>& b, vector<double>& x, int & iter, double & error);
    
};

#endif // COEFFICIENT_MATRIX_HPP