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

#include "../include/GradientGreenGauss.hpp"

void GradientGreenGauss::computeGradient (SolutionSnapshot *solutionSnapshot) {
  
  #pragma omp parallel for default(shared)
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    
    // Declare private variables for OpenMP
    double ownerElementFaceInverseWeight, elementVolume, *faceDimNormal, *faceNonDimNormal;
    double tmpVelocityComponent;
    
    Element *ownerElement, *neighborElement;
    Face *currentFace;
    
    double tangentialVelocity[DIM_BUFFER];
    
    // Main loop
    ownerElement = mesh_->getElement(iElement);
    
    // Clear solution
    for (short iVariable=0; iVariable<nVariable_; iVariable++) {
      for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
        variableGradient_[iElement][iVariable][iDim] = 0;
      }
    }
    
    // Calculate fluxes
    for (int iFace=0; iFace<ownerElement->getNumFace(); iFace++) {
      currentFace = ownerElement->getSingleFace(iFace);
      faceDimNormal = currentFace->getDimNormal();
      neighborElement = ownerElement->getFaceAdjacentElement(iFace);
      ownerElementFaceInverseWeight = ownerElement->getFaceInverseDistanceWeight(iFace);
      
      for (short iVariable=0; iVariable<nVariable_; iVariable++) {        
        for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
          
          // External face (assume value on face same as cell center)
          if (neighborElement == NULL) {
            variableGradient_[iElement][iVariable][iDim] += faceDimNormal[iDim]*(
              solutionSnapshot->getVariable(iElement, iVariable, variableType_));
          }
          // Internal face (face normal pointing outwards of element)
          else if (currentFace->getOwnerElement()->getId() == iElement) {
            variableGradient_[iElement][iVariable][iDim] += faceDimNormal[iDim]*(
              ownerElementFaceInverseWeight*
                solutionSnapshot->getVariable(iElement, iVariable, variableType_) + 
              (1-ownerElementFaceInverseWeight)*
                solutionSnapshot->getVariable(neighborElement->getId(), iVariable, variableType_)
              );
          }
          // Internal face (face normal pointing inwards of element)
          else {
            variableGradient_[iElement][iVariable][iDim] -= faceDimNormal[iDim]*(
              ownerElementFaceInverseWeight*
                solutionSnapshot->getVariable(iElement, iVariable, variableType_) + 
              (1-ownerElementFaceInverseWeight)*
                solutionSnapshot->getVariable(neighborElement->getId(), iVariable, variableType_)
              );
          }
          
        }
      }
      
      // Correct velocity on wall for slip wall
      if (currentFace->getBoundaryType() == BoundaryType::SLIP_WALL) {
        
        // Compute contravariant velocity
        tmpVelocityComponent = 0;
        faceNonDimNormal = currentFace->getNonDimNormal();
        for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
          tmpVelocityComponent += 
            solutionSnapshot->getVariable(
              iElement, 
              config_->mapVelocityComponent(variableType_, iDim), 
              variableType_)*
            faceNonDimNormal[iDim];
        }
        
        // Compute tangential velocity vector
        for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
          tangentialVelocity[iDim] = 
            solutionSnapshot->getVariable(
              iElement, 
              config_->mapVelocityComponent(variableType_, iDim), 
              variableType_) - 
            tmpVelocityComponent*faceNonDimNormal[iDim];
        }
        
        // Correct previously added flux
        for (short iVariable=0; iVariable<mesh_->getNumDim(); iVariable++) {
          for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
            variableGradient_[iElement][config_->mapVelocityComponent(variableType_, iVariable)][iDim] += (
              tangentialVelocity[iVariable]*faceDimNormal[iDim] - // Subtract previously added flux
              solutionSnapshot->getVariable(
                iElement, 
                config_->mapVelocityComponent(variableType_, iVariable), 
                variableType_)*
              faceDimNormal[iDim]
            );
          }
        }
        
      }
      
    } // for iFace (in each element)
    
    // Calculate gradient
    elementVolume = ownerElement->getVolume();
    for (short iVariable=0; iVariable<nVariable_; iVariable++) {
      for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
        variableGradient_[iElement][iVariable][iDim] /= elementVolume;
      }
    }
    
  } // #pragma omp parallel for iElement
  
}