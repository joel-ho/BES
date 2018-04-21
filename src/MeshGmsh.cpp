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

#include "../include/MeshGmsh.hpp"

// Gmsh element 1 as face
void FaceGmshOne::initialize () {
  initializeBase_(2, 2);
}

void FaceGmshOne::computeCenter () {
  center_->initialize(nDim_);
  for (int i=0; i<nDim_; i++) {
    center_->setSingleCoord(
    0.5*(vertices_[0]->getSingleCoord(i)+vertices_[1]->getSingleCoord(i)), 
    i);
  }
}

void FaceGmshOne::computeNormal () {
  dimNormal_[0] = vertices_[1]->getSingleCoord(1) - vertices_[0]->getSingleCoord(1);
  dimNormal_[1] = vertices_[0]->getSingleCoord(0) - vertices_[1]->getSingleCoord(0);
  area_ = MeshMath::computeMagnitude(dimNormal_, nDim_);
  for (short i=0; i<nDim_; i++) {
    nonDimNormal_[i] = dimNormal_[i]/area_;
  }
}

// Gmsh element 2 as face
void FaceGmshTwo::initialize () {
  initializeBase_(3, 3);
}

void FaceGmshTwo::computeCenter () {
  center_->initialize(nDim_);
  for (short i=0; i<nDim_; i++) {
    center_->setSingleCoord(
      (vertices_[0]->getSingleCoord(i) +
      vertices_[1]->getSingleCoord(i) + 
      vertices_[2]->getSingleCoord(i))/3, i);
  }
}

void FaceGmshTwo::computeNormal () {
  MeshMath::crossProduct(dimNormal_, vertices_[0], vertices_[1], vertices_[2]);
  for (short i=0; i<nDim_; i++) {
    dimNormal_[i] /= 2.0;
  }
  area_ = MeshMath::computeMagnitude(dimNormal_, nDim_);
  for (short i=0; i<nDim_; i++) {
    nonDimNormal_[i] = dimNormal_[i]/area_;
  }
}

// Gmsh element 2 as element
void ElementGmshTwo::initialize () {
  initializeBase_(3, 3, 2);
  vtkType_ = 5;
}

void ElementGmshTwo::computeVolume () {
  volume_ = 0.5*abs(
    (vertices_[0]->getSingleCoord(0) - vertices_[1]->getSingleCoord(0))*
    (vertices_[2]->getSingleCoord(1) - vertices_[1]->getSingleCoord(1)) - 
    (vertices_[2]->getSingleCoord(0) - vertices_[1]->getSingleCoord(0))*
    (vertices_[0]->getSingleCoord(1) - vertices_[1]->getSingleCoord(1)) );
}

void ElementGmshTwo::computeCentroid () {
  centroid_ = new Vertex;
  centroid_->initialize(nDim_);
  for (short i=0; i<nDim_; i++) {
    centroid_->setSingleCoord(
      (vertices_[0]->getSingleCoord(i) +
      vertices_[1]->getSingleCoord(i) + 
      vertices_[2]->getSingleCoord(i))/3, i);
  }
}

// Gmsh element 4 as element
void ElementGmshFour::initialize () {
  initializeBase_(4, 4, 3);
  vtkType_ = 10;
}

void ElementGmshFour::computeVolume () {
  
  short tmpVertexOrder[] = {
    0, 2, 1, 
    1, 2, 3, 
    0, 3, 2, 
    0, 1, 3};
  double tmpFaceCenter[3], tmpDimNormal[3];
  
  volume_ = 0;
  for (short iFace=0; iFace<4; iFace++) {
    
    MeshMath::crossProduct(tmpDimNormal, 
      vertices_[tmpVertexOrder[iFace*3+0]], 
      vertices_[tmpVertexOrder[iFace*3+1]], 
      vertices_[tmpVertexOrder[iFace*3+2]]);
    
    for (short iDim=0; iDim<3; iDim++) {
      tmpDimNormal[iDim] /= 2.0;
      
      tmpFaceCenter[iDim] = 0;
      for (short iVertex=0; iVertex<3; iVertex++) {
        tmpFaceCenter[iDim] += vertices_[ tmpVertexOrder[iFace*3+iVertex] ]->getSingleCoord(iDim);
      }
      tmpFaceCenter[iDim] /= 3.0;
      
    }
    
    volume_ += MeshMath::dotProduct(tmpFaceCenter, tmpDimNormal, 3);
    
  }

}

void ElementGmshFour::computeCentroid () {
  centroid_ = new Vertex;
  centroid_->initialize(nDim_);
  for (short i=0; i<nDim_; i++) {
    centroid_->setSingleCoord(
      (vertices_[0]->getSingleCoord(i) +
      vertices_[1]->getSingleCoord(i) + 
      vertices_[2]->getSingleCoord(i) + 
      vertices_[2]->getSingleCoord(i))/4, i);
  }
}