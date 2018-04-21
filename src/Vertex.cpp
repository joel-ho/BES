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

#include "../include/Vertex.hpp"

Vertex::Vertex () {
  coords_ = NULL;
}

Vertex::~Vertex() {
  if (coords_ != NULL) {
    delete [] coords_;
    coords_ = NULL;
  }
}

void Vertex::initialize (double *coords, const short &nDim) {
  coords_ = new double [nDim];
  for (short i = 0; i< nDim ; i++) {
    coords_[i] = coords[i];
  }
}

void Vertex::initialize (const short &nDim) {
  coords_ = new double [nDim];
}

void Vertex::setId(const unsigned long &id) {
  id_ = id;
}

void Vertex::setCoords (double *coords, const short &nDim) {
  for (short i = 0; i< nDim ; i++) {
    coords_[i] = coords[i];
  }
}

void Vertex::setSingleCoord (double singleCoord, const short &iCoord) {
  coords_[iCoord] = singleCoord;
}

unsigned long Vertex::getId() {
  return id_;
}

void Vertex::getCoords(double *resultantCoords, const short &nDim) {
  for (short i = 0; i< nDim ; i++) {
    resultantCoords[i] = coords_[i];
  }
}

void Vertex::getSingleCoord(double &singleCoords, const short &iCoord) {
  singleCoords = coords_[iCoord];
}

double Vertex::getSingleCoord(const short &iCoord) {
  return coords_[iCoord];
}

void Vertex::printCoords(const short &nDim) {
  cout << "[";
  for (short i=0; i<nDim-1; i++) {
    cout << coords_[i] << ", ";
  }
  cout << coords_[nDim-1] << "]";
}
