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

#include "../include/Solver.hpp"

Solver::Solver () {}
Solver::~Solver () {
  delete solution_;
}

void Solver::initializeBase_ (Config *config, Gas *gas, Mesh *mesh) {
  
  config_ = config;
  gas_ = gas;
  mesh_ = mesh;
  
  solution_ = new Solution;
  solution_->initialize(config_, gas_, mesh_);
  
  gradient_ = NULL;
  limiter_ = NULL;
  timeStepCalculator_ = NULL;
  
}

void Solver::getFirstOrderFaceVariable_(
  SolutionSnapshot *solutionSnapshot,
  unsigned long &iFace,
  double * leftVariable, 
  double * rightVariable) 
{
  
  Element *element;
  
  element = mesh_->getFace(iFace)->getOwnerElement();
  for (short iVar=0; iVar<config_->getNumConservedVariable(); iVar++) {
    leftVariable[iVar] = solutionSnapshot->getConservedVariable(element->getId(), iVar);
  }
  
  element = mesh_->getFace(iFace)->getNeighborElement();
  if (element != NULL) { 
    for (short iVar=0; iVar<config_->getNumConservedVariable(); iVar++) {
      rightVariable[iVar] = solutionSnapshot->getConservedVariable(element->getId(), iVar);
    }
  }
  else {
    for (short iVar=0; iVar<config_->getNumConservedVariable(); iVar++) {
      rightVariable[iVar] = 0;
    }
  }
  
}

void Solver::computeFaceConservedVariable_ (
  SolutionSnapshot *solutionSnapshot,
  Gradient *gradient,
  Limiter *limiter,
  unsigned long &iFace,
  double *leftVariable, 
  double *rightVariable) {
  
  double tmpLimiterValue;
  Element *tmpElement;
  
  tmpElement = mesh_->getFace(iFace)->getOwnerElement();
  if (tmpElement != NULL) {
    for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
      tmpLimiterValue = 
        limiter->getLimiterValue(tmpElement->getId(), iVariable);
      leftVariable[iVariable] = 
        solutionSnapshot->getConservedVariable(tmpElement->getId(), iVariable);
      for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
        leftVariable[iVariable] += 
          tmpLimiterValue*
          gradient->getGradient(tmpElement->getId(), iVariable, iDim)*
          tmpElement->getVectorToFaceFromFaceId(iFace, iDim);
      }
    }
  }
  else {
    cout << "ERROR: Encountered face with no owner element. Exiting..." << endl;
    exit(EXIT_FAILURE);
  }
  
  tmpElement = mesh_->getFace(iFace)->getNeighborElement();
  if (tmpElement != NULL) {
    for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
      tmpLimiterValue = 
        limiter->getLimiterValue(tmpElement->getId(), iVariable);
      rightVariable[iVariable] = 
        solutionSnapshot->getConservedVariable(tmpElement->getId(), iVariable);
      for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
        rightVariable[iVariable] += 
          tmpLimiterValue*
          gradient->getGradient(tmpElement->getId(), iVariable, iDim)*
          tmpElement->getVectorToFaceFromFaceId(iFace, iDim);
      }
    }
  }
  else { // External boundary
    for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
      rightVariable[iVariable] = 0;
    }
  }
  
}

void Solver::writeToFile_(const unsigned long &iTimeStep, const short &iCurrentSnapshot) {
  solution_->writeSolution(iTimeStep, iCurrentSnapshot);
  if (config_->getOutputSlipWall() == OutputSlipWall::YES) {
    solution_->writeSlipWall(iTimeStep, iCurrentSnapshot);
  }
}

void Solver::appendToResidual_(string residualFileName, string residual) {
  
  ofstream residualFile;
  
  residualFile.open(residualFileName, ios::app);
  residualFile << residual;
  residualFile.close();
  
}