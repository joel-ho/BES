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

#include "../include/MatrixSparse.hpp"

MatrixSparse::MatrixSparse () {}

MatrixSparse::~MatrixSparse () {}

void MatrixSparse::initialize (ptrdiff_t nRow, ptrdiff_t nNonZero) {
  nRow_ = nRow;
  currentRow_ = 0;
  currentNonZero_ = 0;
  row_.resize(nRow+1, 0);
  col_.resize(nNonZero, 0);
  val_.resize(nNonZero, 0);
}

void MatrixSparse::insertVal (ptrdiff_t iCol, double val) {
  col_[currentNonZero_] = iCol;
  val_[currentNonZero_] = val;
  currentNonZero_++;
}

void MatrixSparse::closeRow () {
  row_[currentRow_+1] = currentNonZero_;
  currentRow_++;
}

ptrdiff_t MatrixSparse::getIdx_ (ptrdiff_t iRow, ptrdiff_t iCol) {
  // Returns -1 if entry not in matrix
  
  vector<ptrdiff_t>::iterator colIt = 
    find(col_.begin()+row_[iRow], col_.begin()+row_[iRow+1], iCol);
  
  // Vector iterator satisfies RandomAccessIterator requirements. Meaning
  // minus operator returns a vector<type>::iterator::difference_type object
  // which is ptrdiff_t in this case.
  return ( colIt != col_.begin()+row_[iRow+1] )? colIt-col_.begin():-1;
  
}

double MatrixSparse::getVal(ptrdiff_t iRow, ptrdiff_t iCol) {
  ptrdiff_t idx;
  idx = getIdx_(iRow, iCol);
  return (idx != -1)? val_[idx]:0;
}

void MatrixSparse::changeVal(ptrdiff_t iRow, ptrdiff_t iCol, double newVal) {
  ptrdiff_t idx;
  idx = getIdx_(iRow, iCol);
  if (idx != -1) {
    val_[idx] = newVal;
  }
}

void MatrixSparse::addToVal(ptrdiff_t iRow, ptrdiff_t iCol, double addVal) {
  ptrdiff_t idx;
  idx = getIdx_(iRow, iCol);
  if (idx != -1) {
    #pragma omp atomic
    val_[idx] += addVal;
  }
}

void MatrixSparse::changeAllVal(double newVal) {
  for (
    vector<double>::iterator it=val_.begin(); 
    it!=val_.end(); 
    ++it
  ){
    *it = newVal;
  }
}

ptrdiff_t MatrixSparse::getNumRow() {
  return nRow_;
}

ptrdiff_t MatrixSparse::getNumNonZero() {
  return val_.size();
}

ptrdiff_t * MatrixSparse::getRowPtr () {
  return row_.data();
}

ptrdiff_t * MatrixSparse::getColPtr () {
  return col_.data();
}

double * MatrixSparse::getValPtr () {
  return val_.data();
}

void MatrixSparse::printSquareMatrix() {
  for (ptrdiff_t iRow=0; iRow<nRow_; iRow++) {
    for (ptrdiff_t iCol=0; iCol<nRow_; iCol++) {
      cout << getVal(iRow, iCol) << ", ";
    }
    cout << endl;
  }
}