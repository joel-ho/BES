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

#pragma once
#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <iostream>
#include <map>

#include "../include/Options.hpp"
#include "../include/Mesh.hpp"

using namespace std;

class Vertex;
class Face;

class Element {
  
  protected:
    
    short nDim_, vtkType_;
    int nVertex_, nFace_;
    unsigned long id_;
    double volume_;
    double **vectorToFace_;
    double *faceInverseDistanceWeight_;
    
    map <unsigned long, int> faceGlobalToLocalMap_;
    
    Vertex *centroid_;
    
    // Pointers to external entities in the mesh 
    // (do not call delete on object pointed to)
    Vertex **vertices_;
    Face **faces_;
    Element **faceAdjacentElement_;
    
    void initializeBase_(int nVertex, int nFace, const short &nDim);
    
    void printVectorToFace_(const short &iFace);
    
  public:
    
    Element();
    virtual ~Element();
    
    virtual void initialize() = 0;
    
    virtual void computeCentroid() = 0;
    virtual void computeVolume() = 0;
    
    void setId(unsigned long &id);
    void setSingleVertex(Vertex *vertex, const int &iVertex);
    void setSingleFace(Face *face, const int &iFace);
    void setFaceAdjacentElement(Element* element, const int& iFace);
    
    void computeAllVectorToFace();
    void computeAllFaceInverseDistanceWeightage();
    
    double getVolume();
    void getVolume(double &elementVolume);
    
    Vertex* getCentroid();
    void getCentroid(Vertex *elementCentroid);
    
    Vertex* getSingleVertex(const int &iVertex);
    void getSingleVertex(Vertex *elementVertex, const int &iVertex);
    
    Face* getSingleFace(const int &iFace);
    
    Element* getFaceAdjacentElement(const int& iFace);
    
    double* getVectorToFace(const int &iFace);
    double* getVectorToFaceFromFaceId(const unsigned long &faceId);
    double getVectorToFaceFromFaceId(const unsigned long &faceId, const short &iDim);
    
    double getFaceInverseDistanceWeight(const int& iFace);
    
    int getNumVertex();
    int getNumFace();
    
    unsigned long getId();
    
    short getVtkType();
    
    void printProperties();
    
};

#endif // ELEMENT_HPP
