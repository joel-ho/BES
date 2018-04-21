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

#include "../include/Element.hpp"

Element::Element () {
  
  vectorToFace_ = NULL;
  faceInverseDistanceWeight_ = NULL;
  centroid_ = NULL;
  
  faceGlobalToLocalMap_.clear();
  
  vertices_ = NULL;
  faces_ = NULL;
  faceAdjacentElement_ = NULL;
  
}

Element::~Element () {
  
  if (vectorToFace_ != NULL) { 
    for (int iFace=0; iFace<nFace_; iFace++) {
      delete [] vectorToFace_[iFace];
    }
    delete [] vectorToFace_;
    vectorToFace_ = NULL;
  }
  
  if (faceInverseDistanceWeight_ != NULL) { 
    delete [] faceInverseDistanceWeight_;
    faceInverseDistanceWeight_ = NULL;
  }
  
  if (centroid_ != NULL) { 
    delete centroid_;
    centroid_ = NULL;
  }
  
  if (vertices_ != NULL) { 
    delete [] vertices_;
    vertices_ = NULL;
  }
  
  if (faces_ != NULL) { 
    delete [] faces_;
    faces_ = NULL;
  }
  
  if (faceAdjacentElement_ != NULL) { 
    delete [] faceAdjacentElement_;
    faceAdjacentElement_ = NULL;
  }
  
}

void Element::initializeBase_ (int nVertex, int nFace, const short &nDim) {
  
  nDim_ = nDim;
  nVertex_ = nVertex;
  nFace_ = nFace;
  
  vertices_ = new Vertex* [nVertex_];
  for (int iVertex=0; iVertex<nVertex_; iVertex++) {
    vertices_[iVertex] = NULL;
  }
  
  faces_ = new Face* [nFace_];
  faceAdjacentElement_ = new Element* [nFace_];
  faceInverseDistanceWeight_ = new double [nFace_];
  for (int iFace=0; iFace<nFace_; iFace++) {
    faces_[iFace] = NULL;
    faceAdjacentElement_[iFace] = NULL;
    faceInverseDistanceWeight_[iFace] = 0;
  }
  
}

void Element::setId(unsigned long &id) {
  id_ = id;
}

void Element::setSingleVertex(Vertex *vertex, const int &iVertex) {
  vertices_[iVertex] = vertex;
}

void Element::setSingleFace(Face *face, const int &iFace) {
  faces_[iFace] = face;
  faceGlobalToLocalMap_[face->getId()] = iFace;
}

void Element::setFaceAdjacentElement(Element* element, const int& iFace) {
  faceAdjacentElement_[iFace] = element;
}

void Element::computeAllVectorToFace() {
  vectorToFace_ = new double* [nFace_];
  for (int iFace=0; iFace<nFace_; iFace++) {
    vectorToFace_[iFace] = new double [nDim_];    
    MeshMath::computeVector(centroid_, faces_[iFace]->getCenter(), 
      vectorToFace_[iFace], nDim_);
  }
}

void Element::computeAllFaceInverseDistanceWeightage () {
  
  double ownerFaceDistance, neighborFaceDistance;
  
  for (int iFace=0; iFace<nFace_; iFace++) {
    if (faceAdjacentElement_[iFace] != NULL) {
      ownerFaceDistance = MeshMath::computeDistance(
        centroid_,
        faces_[iFace]->getCenter(),
        nDim_);
      neighborFaceDistance = MeshMath::computeDistance(
        faceAdjacentElement_[iFace]->getCentroid(),
        faces_[iFace]->getCenter(),
        nDim_);
        
      faceInverseDistanceWeight_[iFace] = neighborFaceDistance/
        (ownerFaceDistance+neighborFaceDistance);
        
    }
  }
}


double Element::getVolume() {
  return volume_;
}

void Element::getVolume (double &elementVolume) {
  elementVolume = volume_;
}

Vertex* Element::getCentroid() {
  return centroid_;
}

void Element::getCentroid (Vertex *elementCentroid) {
  elementCentroid = centroid_;
}

Vertex* Element::getSingleVertex (const int &iVertex) {
  Vertex *result;
  result = vertices_[iVertex];
  return result;
}

void Element::getSingleVertex (Vertex *elementVertex, const int &iVertex) {
  elementVertex = vertices_[iVertex];
}

Face* Element::getSingleFace(const int& iFace) {
  return faces_[iFace];
}

Element* Element::getFaceAdjacentElement(const int& iFace) {
  return faceAdjacentElement_[iFace];
}

double* Element::getVectorToFace(const int &iFace) {
  return vectorToFace_[iFace];
}

double* Element::getVectorToFaceFromFaceId(const unsigned long &faceId) {
  return vectorToFace_[faceGlobalToLocalMap_[faceId]];
}

double Element::getVectorToFaceFromFaceId(const unsigned long &faceId, const short &iDim) {
  return vectorToFace_[faceGlobalToLocalMap_[faceId]][iDim];
}

double Element::getFaceInverseDistanceWeight (const int &iFace) {
  return faceInverseDistanceWeight_[iFace];
}

int Element::getNumVertex() {
  return nVertex_;
}

int Element::getNumFace() {
  return nFace_;
}

unsigned long Element::getId() {
  return id_;
}

void Element::printVectorToFace_(const short &iFace) {
  cout << "[";
  for (short i=0; i<nDim_-1; i++) {
    cout << vectorToFace_[iFace][i] << ", ";
  }
  cout << vectorToFace_[iFace][nDim_-1] << "]";
}

short Element::getVtkType () {
  return vtkType_;
}

void Element::printProperties() {
  
  cout << "Element ID: " << id_ << endl;
  cout << "  Number of vertices: " << nVertex_ << ", Number of faces: " << nFace_ << endl;
  cout << "  Volume: " << volume_ << endl;
  cout << "  Centroid coordinates: ";
  centroid_->printCoords(nDim_);
  cout << endl;
  cout << "  Vertices coordinates: " << endl;
  for (short iVertex=0; iVertex<nVertex_; iVertex++) {
    cout << "    ";
    vertices_[iVertex]->printCoords(nDim_);
    cout << endl;
  }
  cout << "  Face nondimensional normal vectors: " << endl;
  for (short iFace=0; iFace<nFace_; iFace++) {
    cout << "    Face owner element: " << faces_[iFace]->getOwnerElement()->getId() << ", vector: ";
    faces_[iFace]->printNonDimNormal();
    cout << endl;
  }
  cout << "  Face dimensional normal vectors: " << endl;
  for (short iFace=0; iFace<nFace_; iFace++) {
    cout << "    Face owner element: " << faces_[iFace]->getOwnerElement()->getId() << ", vector: ";
    faces_[iFace]->printDimNormal();
    cout << endl;
  }
  cout << "  Inverse distance weights: " << endl;
  for (short iFace=0; iFace<nFace_; iFace++) {
    cout << "    Face owner element: " << faces_[iFace]->getOwnerElement()->getId() << ", weight: ";
    cout << faceInverseDistanceWeight_[iFace];
    cout << endl;
  }
  cout << "  Vector to face: " << endl;
  for (short iFace=0; iFace<nFace_; iFace++) {
    cout << "    Face: " << iFace << ", vector: ";
    printVectorToFace_(iFace);
    cout << endl;
  }
  cout << "  Face adjacent elements: " << endl;
  for (short iFace=0; iFace<nFace_; iFace++) {
    if (faceAdjacentElement_[iFace]!=NULL) {
      cout << "    Face: " << iFace << ", element: " << faceAdjacentElement_[iFace]->getId();
    }
    else {
      cout << "    Face: " << iFace << ", no adjacent element.";
    }
    cout << endl;
  }
}
