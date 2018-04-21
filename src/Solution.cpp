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

#include "../include/Solution.hpp"

Solution::Solution() {
  
  nCurrentSnapshots_ = 0;
  
  currentSnapshots_ = NULL;
  
  config_ = NULL;
  gas_ = NULL;
  mesh_ = NULL;
  
}

Solution::~Solution() {
  if (currentSnapshots_ != NULL) {
    for (short iSol=0; iSol<nCurrentSnapshots_; iSol++) {
      delete currentSnapshots_[iSol];
    }
    delete [] currentSnapshots_;
    currentSnapshots_ = NULL;
  }
}

void Solution::initialize(Config* config, Gas *gas, Mesh* mesh) {
  
  // Save class properties
  config_ = config;
  gas_ = gas;
  mesh_ = mesh;
  
  // Set some properties
  nCurrentSnapshots_ = config->getNumCurrentSolution();
  
  // Set up storage
  currentSnapshots_ = new SolutionSnapshot* [nCurrentSnapshots_];
  
  for (short iSol=0; iSol<nCurrentSnapshots_; iSol++) {
    currentSnapshots_[iSol] = new SolutionSnapshot();
    currentSnapshots_[iSol]->initialize(config_, mesh_, gas_);
  }
  
}

inline short Solution::indexCurrentSnapshot_ (
  const unsigned long &iTimeStep, const short &iCurrentSnapshot) {
  
  return (
    ((iTimeStep%2)*(nCurrentSnapshots_-1))%nCurrentSnapshots_ + 
    ((iTimeStep%2)*(-2) + 1)*iCurrentSnapshot
  );
  
}

void Solution::setConservedVariable(
  const unsigned long &iTimeStep, const short &iCurrentSnapshot,
  const unsigned long &iElement, double *conservedVariable) {
  currentSnapshots_[indexCurrentSnapshot_(iTimeStep, iCurrentSnapshot)]->
    setConservedVariable(iElement, conservedVariable);
  
}

void Solution::setConservedVariable(
  const unsigned long &iTimeStep, const short &iCurrentSnapshot,
  const unsigned long &iElement, const short &iVariable,
  double value) {
    
  currentSnapshots_[indexCurrentSnapshot_(iTimeStep, iCurrentSnapshot)]->
    setConservedVariable(iElement, iVariable, value);
  
}

double* Solution::getConservedVariable(
  const unsigned long &iTimeStep, const short &iCurrentSnapshot, 
  const unsigned long &iElement) {
  
  return (
    currentSnapshots_[indexCurrentSnapshot_(iTimeStep, iCurrentSnapshot)]->
      getConservedVariable(iElement)
  );
  
}

double Solution::getConservedVariable(
  const unsigned long &iTimeStep, const short &iCurrentSnapshot, 
  const unsigned long &iElement, const short &iVariable) {
  
  return (
    currentSnapshots_[indexCurrentSnapshot_(iTimeStep, iCurrentSnapshot)]->
      getConservedVariable(iElement, iVariable)
  );
  
}

double Solution::getPrimitiveVariable(
  const unsigned long &iTimeStep, const short &iCurrentSnapshot, 
  const unsigned long &iElement, const short &iVariable) {
  
  return (
    currentSnapshots_[indexCurrentSnapshot_(iTimeStep, iCurrentSnapshot)]->
      getPrimitiveVariable(iElement, iVariable)
  );
  
}

SolutionSnapshot* Solution::getCurrentSnapshot(
  const unsigned long &iTimeStep, const short &iCurrentSnapshot) {
    return currentSnapshots_[indexCurrentSnapshot_(iTimeStep, iCurrentSnapshot)];
}

void Solution::writeSolution (const unsigned long &iTimeStep, const short &iCurrentSnapshot) {
  
  ostringstream outputNameStream;
  
  /*Solution->*/getCurrentSnapshot(iTimeStep, iCurrentSnapshot)->computeAllPrimitiveVariable();
  
  outputNameStream << config_->getOutputFolderPath() << "/solution_" << iTimeStep+1;
  switch (config_->getOutputFormat()) {
  
    case OutputFormat::VTK:
      outputNameStream << ".vtk";
      getCurrentSnapshot(iTimeStep, iCurrentSnapshot)->writePrimitiveToVtk(outputNameStream.str());
      break;
      
    case OutputFormat::CSV:
      outputNameStream << ".csv";
      getCurrentSnapshot(iTimeStep, iCurrentSnapshot)->writePrimitiveToCsv(outputNameStream.str());
      break;
      
    default:
      break;
  
  }
  
}

void Solution::writeSlipWall (const unsigned long &iTimeStep, const short &iCurrentSnapshot) {
  
  ostringstream outputNameStream;
  
  outputNameStream << config_->getOutputFolderPath() << "/slip_wall_" << iTimeStep+1 << ".csv";
  /*Solution->*/getCurrentSnapshot(iTimeStep, iCurrentSnapshot)->writeSlipWallToCsv(outputNameStream.str());
  
}