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

#include "../include/TimeStepCalculatorLocal.hpp"

void TimeStepCalculatorLocal::initialize (
  Config *config, Gas *gas, Mesh *mesh) {
  
  initializeBase_(config, gas, mesh);
  
  allTimeStep_ = new double [mesh_->getNumElement()];
  projectedLength_ = new double* [mesh_->getNumElement()];
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    
    allTimeStep_[iElement] = 0;
    
    projectedLength_[iElement] = new double [mesh_->getNumDim()];
    for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
      projectedLength_[iElement][iDim] = 0;
    }
    
  }
  
  computeProjectedLength_();
  
}

TimeStepCalculatorLocal::~TimeStepCalculatorLocal () {
  delete [] allTimeStep_;
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    delete [] projectedLength_[iElement];
  }
  delete [] projectedLength_;
}

void TimeStepCalculatorLocal::computeProjectedLength_() {
  
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    for (int iFace=0; iFace<mesh_->getElement(iElement)->getNumFace(); iFace++) {
      for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
        projectedLength_[iElement][iDim] +=  abs(
          mesh_->getElement(iElement)->getSingleFace(iFace)->getDimNormal(iDim));
      }
    }
    for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
      projectedLength_[iElement][iDim] /= 2.0;
    }
  }
  
}

void TimeStepCalculatorLocal::computeAllTimeStep (
  SolutionSnapshot *solutionSnapshot, const double &courantNum) {
  
  double spectralRadii, soundSpeed, velocity;
  
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    
    spectralRadii = 0.0;
    soundSpeed = gas_->computeSoundSpeed(
      solutionSnapshot->getConservedVariable(iElement));
    
    for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
      
      velocity = abs(
        gas_->computeVelocity(
          solutionSnapshot->getConservedVariable(iElement), iDim)
      );
      
      spectralRadii += (soundSpeed + velocity)*projectedLength_[iElement][iDim];
      
    }
    
    allTimeStep_[iElement] = 
      courantNum*mesh_->getElement(iElement)->getVolume()/spectralRadii;
    
  }
}

double TimeStepCalculatorLocal::getTimeStep(const unsigned long &iElement) {
  return allTimeStep_[iElement];
}