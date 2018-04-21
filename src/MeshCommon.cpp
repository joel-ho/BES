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

#include "../include/MeshCommon.hpp"

double MeshMath::computeDistance (
  double *coordsOne, double *coordsTwo, const short &nDim) {
  
  double distance=0;
  for (short i=0; i<nDim; i++) {
    distance += 
      pow(coordsOne[i] - coordsTwo[i], 2.0);
  }
  distance = sqrt(distance);
  return distance;
}

double MeshMath::computeDistance (
  Vertex *vertexOne, Vertex *vertexTwo, const short &nDim) {
  
  double distance=0;
  for (short i=0; i<nDim; i++) {
    distance += 
      pow(vertexOne->getSingleCoord(i) - vertexTwo->getSingleCoord(i), 2.0);
  }
  distance = sqrt(distance);
  return distance;
}

double MeshMath::computeMagnitude (double *vector, const short &nDim) {
  double magnitude=0;
  for (short i=0; i<nDim; i++) {
    magnitude += pow(vector[i], 2.0);
  }
  magnitude = sqrt(magnitude);
  return magnitude;
}

void MeshMath::nonDimVector (double *vector, const short &nDim, double area) {
  for (short i=0; i<nDim; i++) {
    vector[i] = vector[i]/area;
  }
}

void MeshMath::computeTwoDimOutwardNormal (
  double *result, Vertex *vertexOne, Vertex *vertexTwo) {
  // Returns normal pointing outward of 2D element with vertex number increasing
  // in counter clockwise direction. vertexTwo is in clockwise direction 
  // relative to vertexOne.
  result[0] = vertexTwo->getSingleCoord(1) - vertexOne->getSingleCoord(1);
  result[1] = - vertexTwo->getSingleCoord(0) + vertexOne->getSingleCoord(0);
}

double MeshMath::dotProduct (double *vecOne, double* vecTwo, const short &nDim) {
  double results=0;
  for (short i=0; i<nDim; i++) {
    results += vecOne[i]*vecTwo[i];
  }
  return results;
}

void MeshMath::crossProduct (
  double *result, Vertex *vertexOne, Vertex *vertexTwo, Vertex *vertexThree) {
  // cross v_{3, 2} with v_{1, 2}
  // vector v_{b, a} points from a to b
  result[0] = 
    (vertexThree->getSingleCoord(1) - vertexTwo->getSingleCoord(1))*
    (vertexOne->getSingleCoord(2) - vertexTwo->getSingleCoord(2)) -
    (vertexOne->getSingleCoord(1) - vertexTwo->getSingleCoord(1))*
    (vertexThree->getSingleCoord(2) - vertexTwo->getSingleCoord(2));
  result[1] = 
    (vertexOne->getSingleCoord(0) - vertexTwo->getSingleCoord(0))*
    (vertexThree->getSingleCoord(2) - vertexTwo->getSingleCoord(2)) -
    (vertexThree->getSingleCoord(0) - vertexTwo->getSingleCoord(0))*
    (vertexOne->getSingleCoord(2) - vertexTwo->getSingleCoord(2));
  result[2] = 
    (vertexThree->getSingleCoord(0) - vertexTwo->getSingleCoord(0))*
    (vertexOne->getSingleCoord(1) - vertexTwo->getSingleCoord(1)) -
    (vertexOne->getSingleCoord(0) - vertexTwo->getSingleCoord(0))*
    (vertexThree->getSingleCoord(1) - vertexTwo->getSingleCoord(1));
}

void MeshMath::computeVector(Vertex* vertexOne, Vertex* vertexTwo, 
  double* resultantVector, const short &nDim) {
  // Resultant vector points from vertex one to vertex two
  for (short iDim=0; iDim<nDim; iDim++) {
    resultantVector[iDim] = 
      vertexTwo->getSingleCoord(iDim) - vertexOne->getSingleCoord(iDim);
  }
}
