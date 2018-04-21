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
#ifndef MATRIX_SPARSE_HPP
#define MATRIX_SPARSE_HPP

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class MatrixSparse {
  
  private:
    
    // ptrdiff_t may not be able to count up to max unsigned long 
    // (might run into compatibility issues with unsigned long)
    ptrdiff_t nRow_, currentRow_, currentNonZero_;
    vector<ptrdiff_t> row_, col_; 
    vector<double> val_;
    
    ptrdiff_t getIdx_(ptrdiff_t iRow, ptrdiff_t iCol);
    
  public:
    
    MatrixSparse();
    ~MatrixSparse();
    
    void initialize(ptrdiff_t nRow, ptrdiff_t nNonZero);
    
    // Values can only be added in increasing row order.
    // Once closeRow is called, values will be added to new row.
    void insertVal(ptrdiff_t iCol, double val);
    void closeRow();
    
    // No bounds checking
    double getVal(ptrdiff_t iRow, ptrdiff_t iCol);
    
    // Operations are ignored if value does not exist in iRow, iCol.
    // No bounds checking
    void changeVal(ptrdiff_t iRow, ptrdiff_t iCol, double newVal);
    void addToVal(ptrdiff_t iRow, ptrdiff_t iCol, double addVal);
    void changeAllVal(double newVal);
    
    // Return vectors for linear solver.
    // Pointer to first entry of vector.
    ptrdiff_t getNumRow();
    ptrdiff_t getNumNonZero();
    ptrdiff_t* getRowPtr();
    ptrdiff_t* getColPtr();
    double* getValPtr();
    
    void printSquareMatrix();
    
};

#endif // MATRIX_SPARSE_HPP