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

#include "../include/Limiter.hpp"

Limiter::Limiter () {
  config_ = NULL;
  mesh_ = NULL;
  gradientLimiter_ = NULL;
}

Limiter::~Limiter () {
  if (gradientLimiter_ != NULL) {
    for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      delete [] gradientLimiter_[iElement];
    }
    delete [] gradientLimiter_;
  }
}

void Limiter::initialize(Config *config, Mesh *mesh, VariableType variableType) {
  
  config_ = config;
  mesh_ = mesh;
  variableType_ = variableType;
  
  switch (variableType_) {
  
    case VariableType::CONSERVED:
      nVariable_ = config_->getNumConservedVariable();
      break;
      
    case VariableType::PRIMITIVE:
      nVariable_ = config_->getNumPrimitiveVariable();
      break;
      
    default:
      break;
  
  }
  
  gradientLimiter_ = new double* [ mesh_->getNumElement() ];
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    gradientLimiter_[iElement] = new double [nVariable_];
    for (short iVar=0; iVar<nVariable_; iVar++) {
      gradientLimiter_[iElement][iVar] = 0;
    }
  }
  
}

double Limiter::getLimiterValue (
  const unsigned long &iElement, const short &iVariable) {
  return gradientLimiter_[iElement][iVariable];
}
