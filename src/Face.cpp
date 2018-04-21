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

#include "../include/Face.hpp"

Face::Face () {
  
  nonDimNormal_ = NULL;
  dimNormal_ = NULL;
  center_ = NULL;
  
  vertices_ = NULL;
  ownerElement_ = NULL;
  neighborElement_ = NULL;
  
}

Face::~Face () {
  if (nonDimNormal_ != NULL) {
    delete [] nonDimNormal_; 
    nonDimNormal_ = NULL;
  }
  if (dimNormal_ != NULL) { 
    delete [] dimNormal_;
    dimNormal_ = NULL;
  }
  if (center_ != NULL) { 
    delete center_;
    center_ = NULL;
  }
  if (vertices_ != NULL) {
    delete [] vertices_;
    vertices_ = NULL;
  }
}

void Face::initializeBase_ (const int &nVertex, const short &nDim) {
  nVertex_ = nVertex;
  nDim_ = nDim;
  vertices_ = new Vertex* [nVertex_];
  center_ = new Vertex;
  nonDimNormal_ = new double [nDim_];
  dimNormal_ = new double [nDim_];
}

void Face::setId(unsigned long id) {
  id_ = id;
}

void Face::setSingleVertex (Vertex *vertex, const int &iVertex) {
  vertices_[iVertex] = vertex;
}

void Face::setVertices (Vertex **vertices) {
  for (int i=0; i<nVertex_; i++) {
    vertices_[i] = vertices[i];
  }
}

void Face::setBoundaryType (BoundaryType boundaryType) {
  boundaryType_ = boundaryType;
}

void Face::setOwnerElement (Element *ownerElement) {
  ownerElement_ = ownerElement;
}

void Face::setNeighborElement (Element *neighborElement) {
  neighborElement_ = neighborElement;
}

void Face::flipNormal () {
  for (short i=0; i<nDim_; i++) {
    nonDimNormal_[i] = -nonDimNormal_[i];
    dimNormal_[i] = -dimNormal_[i];
  }
}

BoundaryType Face::getBoundaryType () {
  return boundaryType_;
}

int Face::queryBoundaryExternal () {
  return ExternalBoundaryMap[boundaryType_];
}

void Face::queryBoundaryExternal (int &boundaryExternal) {
  boundaryExternal = ExternalBoundaryMap[boundaryType_];
}

Element* Face::getOwnerElement() {
  return ownerElement_;
}

void Face::getOwnerElement (Element *ownerElement) {
  ownerElement = ownerElement_;
}

Element* Face::getNeighborElement() {
  return neighborElement_;
}

void Face::getNeighborElement (Element *neighborElement) {
  neighborElement = neighborElement_;
}

double* Face::getDimNormal() {
  return dimNormal_;
}

// Face::getDimNormal (double *normal) {
  // for (short i=0; i<nDim_; i++) {
    // normal[i] = dimNormal_[i];
  // }
// }

double Face::getDimNormal(const short &iDim) {
  return dimNormal_[iDim];
}

double* Face::getNonDimNormal() {
  return nonDimNormal_;
}

double Face::getNonDimNormal(const short &iDim) {
  return nonDimNormal_[iDim];
}

Vertex* Face::getCenter() {
  return center_;
}

double Face::getArea() {
  return area_;
}

void Face::getArea(double &area) {
  area = area_;
}

unsigned long Face::getId() {
  return id_;
}

void Face::printDimNormal() {
  cout << "[";
  for (short i=0; i<nDim_-1; i++) {
    cout << dimNormal_[i] << ", ";
  }
  cout << dimNormal_[nDim_-1] << "]";
}

void Face::printNonDimNormal() {
  cout << "[";
  for (short i=0; i<nDim_-1; i++) {
    cout << nonDimNormal_[i] << ", ";
  }
  cout << nonDimNormal_[nDim_-1] << "]";
}
