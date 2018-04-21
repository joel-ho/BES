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

#include "../include/CoefficientMatrix.hpp"

CoefficientMatrix::CoefficientMatrix () {
  matrix_ = NULL;
  config_ = NULL;
  mesh_ = NULL;
}

CoefficientMatrix::~CoefficientMatrix () {
  if (matrix_ != NULL) {
    delete matrix_;
  }
}

unsigned long CoefficientMatrix::computeIdx_(unsigned long iElement, short iVar) {
  return iElement*blockSize_+iVar;
}

void CoefficientMatrix::printOptions_() {
  cout << "  Initializing matrix system..." << endl;
  
  cout << "    Algebraic multi grid preconditioner coarsener: ";
  switch ( config_->getAmgPrecondCoarsener() ) {
  
    case AmgPrecondCoarsener::AGGREGATION:
      cout << "Aggregation...." << endl;
      break;
      
    case AmgPrecondCoarsener::SMOOTHED_AGGREGATION:
      cout << "Smoothed aggregation..." << endl;
      break;
      
    case AmgPrecondCoarsener::SMOOTHED_AGGREGATION_EMIN:
      cout << "Smoothed aggregation with energy minimization..." << endl;
      break;
      
    case AmgPrecondCoarsener::RUGE_STUBEN:
      cout << "Ruge-Stuben..." << endl;
      break;
      
    default:
      break;
  
  }
  
  cout << "    Algebraic multi grid preconditioner smoother: ";
  switch ( config_->getAmgPrecondSmoother() ) {
  
    case AmgPrecondSmoother::DAMPED_JACOBI:
      cout << "Damped Jacobi..." << endl;
      break;
      
    case AmgPrecondSmoother::ILU0:
      cout << "Incomplete LU with zero fill in..." << endl;
      break;
      
    case AmgPrecondSmoother::SPAI0:
      cout << "Sparse approximate inverse with zero fill in..." << endl;
      break;
      
    case AmgPrecondSmoother::GAUSS_SEIDEL:
      cout << "Gauss-Seidel..." << endl;
      break;
      
    default:
      break;
  
  }
  
  switch ( config_->getLinearSolverType() ) {
  
    case LinearSolver::FGMRES:
      cout << "    Linear solver: FGMRES..." << endl;
      break;
    
    case LinearSolver::GMRES:
      cout << "    Linear solver: GMRES..." << endl;
      break;
    
    case LinearSolver::LGMRES:
      cout << "    Linear solver: LGMRES..." << endl;
      break;
      
    case LinearSolver::BICGSTAB:
      cout << "    Linear solver: BICGSTAB..." << endl;
      break;
    
    default:
      break;
  
  }
  cout << "    Max iterations: " << config_->getLinearSolverNumIter() << "..." << endl;
  cout << "    Abs tolerance: " << config_->getLinearSolverTol() << "... ";
}

void CoefficientMatrix::initialize(Config* config, Mesh* mesh) {
  
  unsigned long nNonZero;
  Element* neighborElement;
  
  config_ = config;
  mesh_ = mesh;
  
  printOptions_();
  
  nBlock_ = mesh_->getNumElement();
  blockSize_ = config_->getNumConservedVariable();
  
  nNonZero = 0;
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    nNonZero++;
    for (int iFace=0; iFace<mesh_->getElement(iElement)->getNumFace(); iFace++) {
      neighborElement = mesh_->getElement(iElement)->getFaceAdjacentElement(iFace);
      if (neighborElement != NULL) {
        nNonZero++;
      }
    }
  }
  nNonZero *= blockSize_*blockSize_;
  
  matrix_ = new MatrixSparse;
  matrix_->initialize(nBlock_*blockSize_, nNonZero);
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    for (short iVarDown=0; iVarDown<blockSize_; iVarDown++) {
      
      for (short iVarAcross=0; iVarAcross<blockSize_; iVarAcross++) {
        matrix_->insertVal(computeIdx_(iElement, iVarAcross), 0);
      }// for iVarAcross
      
      for (int iFace=0; iFace<mesh_->getElement(iElement)->getNumFace(); iFace++) {
        neighborElement = mesh_->getElement(iElement)->getFaceAdjacentElement(iFace);
        if (neighborElement != NULL) {
          for (short iVarAcross=0; iVarAcross<blockSize_; iVarAcross++) {
            matrix_->insertVal(computeIdx_(neighborElement->getId(), iVarAcross), 0);
          } // for iVarAcross
        }
      } // for iFace
      
      matrix_->closeRow();
      
    } // for iVarDown
  } // for iElement
  
  cout << "Done." << endl;
  
}

void CoefficientMatrix::changeAllVal(double newVal) {
  matrix_->changeAllVal(newVal);
}

void CoefficientMatrix::addToDiagonal (double value) {
  for (unsigned long i=0; i<nBlock_*blockSize_; i++) {
    matrix_->addToVal(i, i, value);
  }
}

void CoefficientMatrix::addToVal
(
  unsigned long iElementRow, unsigned long iElementCol, 
  short iVarRow, short iVarCol, 
  double addVal
) 
{
  
  matrix_->addToVal(
    computeIdx_(iElementRow, iVarRow),
    computeIdx_(iElementCol, iVarCol),
    addVal);
}

void CoefficientMatrix::addToBlock
(
  unsigned long iElementRow, unsigned long iElementCol, 
  double multiplier, double * addValRowMajor
)
{
  for (short iRow=0; iRow<blockSize_; iRow++) {
    for (short iCol=0; iCol<blockSize_; iCol++) {
      matrix_->addToVal(
        computeIdx_(iElementRow, iRow),
        computeIdx_(iElementCol, iCol),
        multiplier*addValRowMajor[ iRow*blockSize_ + iCol ]);
    }
  }
}

MatrixSparse* CoefficientMatrix::getMatrixSparse () {
  return matrix_;
}

void CoefficientMatrix::solve(vector<double> & b, vector<double> & x, int & iter, double & error) {
  
  using namespace amgcl;
  
  int tmpIter;
  double tmpError;
  unsigned long nMatrixRow;
  
  boost::property_tree::ptree prm;
  
  // Select solver
  typedef backend::builtin<double> Backend;
  
  typedef make_solver
  <
    runtime::amg<Backend>,
    runtime::iterative_solver<Backend>
  >
  Solver;
  
  // Set parameters
  switch ( config_->getAmgPrecondCoarsener() ) {
  
    case AmgPrecondCoarsener::AGGREGATION:
      prm.put("precond.coarsening.type", runtime::coarsening::aggregation);
      break;
      
    case AmgPrecondCoarsener::SMOOTHED_AGGREGATION:
      prm.put("precond.coarsening.type", runtime::coarsening::smoothed_aggregation);
      break;
      
    case AmgPrecondCoarsener::SMOOTHED_AGGREGATION_EMIN:
      prm.put("precond.coarsening.type", runtime::coarsening::smoothed_aggr_emin);
      break;
      
    case AmgPrecondCoarsener::RUGE_STUBEN:
      prm.put("precond.coarsening.type", runtime::coarsening::ruge_stuben);
      break;
      
    default:
      break;
  
  }
  if (config_->getAmgPrecondCoarsener() != AmgPrecondCoarsener::RUGE_STUBEN) {
    prm.put("precond.coarsening.aggr.block_size", blockSize_);
  }
  
  switch ( config_->getAmgPrecondSmoother() ) {
  
    case AmgPrecondSmoother::DAMPED_JACOBI:
      prm.put("precond.relax.type", runtime::relaxation::damped_jacobi);
      break;
      
    case AmgPrecondSmoother::ILU0:
      prm.put("precond.relax.type", runtime::relaxation::ilu0);
      break;
      
    case AmgPrecondSmoother::SPAI0:
      prm.put("precond.relax.type", runtime::relaxation::spai0);
      break;
      
    case AmgPrecondSmoother::GAUSS_SEIDEL:
      prm.put("precond.relax.type", runtime::relaxation::gauss_seidel);
      break;
      
    default:
      break;
  
  }
  
  switch ( config_->getLinearSolverType() ) {
  
    case LinearSolver::FGMRES:
      prm.put("solver.type", runtime::solver::fgmres);
      break;
    
    case LinearSolver::GMRES:
      prm.put("solver.type", runtime::solver::gmres);
      break;
    
    case LinearSolver::LGMRES:
      prm.put("solver.type", runtime::solver::lgmres);
      break;
      
    case LinearSolver::BICGSTAB:
      prm.put("solver.type", runtime::solver::bicgstab);
      break;
    
    default:
      break;
  
  }
  
  prm.put("solver.maxiter", config_->getLinearSolverNumIter());
  prm.put("solver.tol", config_->getLinearSolverTol());
  
  // Set up solver
  nMatrixRow = matrix_->getNumRow();
  Solver solve( 
    amgcl::adapter::zero_copy(
      nMatrixRow,
      matrix_->getRowPtr(), 
      matrix_->getColPtr(), 
      matrix_->getValPtr()
    ), 
    prm
  );
  
  // Solve system
  boost::tie(tmpIter, tmpError) = solve(b, x);
  
  // Write to iteration and error
  iter = tmpIter;
  error = tmpError;
  
}