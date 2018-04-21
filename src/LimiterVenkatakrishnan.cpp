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

#include "../include/LimiterVenkatakrishnan.hpp"

void LimiterVenkatakrishnan::computeLimiterTerm(
  SolutionSnapshot *solutionSnapshot, Gradient *gradient) {
  
  #pragma omp parallel for default(shared)
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    
    // Declare private variables for OpenMP
    double neighborVariable;
    double deltaOneMin, deltaOneMax, deltaTwo, epsilon, characteristicLength;
    double *vectorToFace;
    
    vector <double> allLimiterValues;
    
    Element *ownerElement, *neighborElement;
    
    allLimiterValues.reserve( 2*mesh_->getNumDim() );
    
    // Main loop
    ownerElement = mesh_->getElement(iElement);
    
    // Calculate epsilon
    characteristicLength = pow(ownerElement->getVolume(), 
      1.0/(mesh_->getNumDim()));
    epsilon = pow(
      (config_->getVenkatakrishnanLimiterConstant()*characteristicLength), 
      (3.0/2.0));
      
    for (short iVariable=0; iVariable<nVariable_; iVariable++) {
      
      deltaOneMin = solutionSnapshot->getVariable(iElement, iVariable, variableType_);
      deltaOneMax = solutionSnapshot->getVariable(iElement, iVariable, variableType_);
      
      // Set results to zero
      gradientLimiter_[iElement][iVariable] = 0;
      allLimiterValues.clear();
      allLimiterValues.push_back(1.0);
      
      // Compute delta_1,min/max
      for (int iFace=0; iFace<ownerElement->getNumFace(); iFace++) {
        neighborElement = ownerElement->getFaceAdjacentElement(iFace);
        if (neighborElement != NULL) {
          neighborVariable = solutionSnapshot->getVariable(neighborElement->getId(), iVariable, variableType_);
          deltaOneMin = (neighborVariable < deltaOneMin)? neighborVariable:deltaOneMin;
          deltaOneMax = (neighborVariable > deltaOneMax)? neighborVariable:deltaOneMax;
        }
      }
      deltaOneMin -= solutionSnapshot->getVariable(iElement, iVariable, variableType_);
      deltaOneMax -= solutionSnapshot->getVariable(iElement, iVariable, variableType_);
      
      // Computer limiter value on each face
      for (int iFace=0; iFace<ownerElement->getNumFace(); iFace++) {
        neighborElement = ownerElement->getFaceAdjacentElement(iFace);
        if (neighborElement != NULL) {
          
          // Compute delta_2
          vectorToFace = ownerElement->getVectorToFace(iFace);
          deltaTwo = 0;
          for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {

            deltaTwo += gradient->getGradient(iElement, iVariable, iDim)*
              vectorToFace[iDim];
          }
          deltaTwo /= 2.0;

          allLimiterValues.push_back( 
            computeFaceLimiter_(deltaOneMin, deltaOneMax, deltaTwo, epsilon));
            
        }
      } // for iFace
      
      // Get minimum limiter term
      gradientLimiter_[iElement][iVariable] = *min_element(
        allLimiterValues.begin(), allLimiterValues.end());
      
    } // for iVariable
  
  } // #pragma omp parallel for iElement
  
}

double LimiterVenkatakrishnan::computeFaceLimiter_(
  double &deltaOneMin, double &deltaOneMax, 
  double &deltaTwo, double &epsilon) {

  if (deltaTwo > 0) {
    return (
      (1/deltaTwo)*
      ((deltaOneMax*deltaOneMax+epsilon*epsilon)*deltaTwo + 
        2.0*deltaTwo*deltaTwo*deltaOneMax)
      /
      (deltaOneMax*deltaOneMax + 2.0*deltaTwo*deltaTwo + 
        deltaOneMax*deltaTwo + epsilon*epsilon)
    );
  }
  else if (deltaTwo < 0) {
    return (
      (1/deltaTwo)*
      ((deltaOneMin*deltaOneMin+epsilon*epsilon)*deltaTwo + 
        2.0*deltaTwo*deltaTwo*deltaOneMin)
      /
      (deltaOneMin*deltaOneMin + 2.0*deltaTwo*deltaTwo + 
        deltaOneMin*deltaTwo + epsilon*epsilon)
    );
  }
  else {
    return 1.0;
  }
}