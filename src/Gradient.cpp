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

#include "../include/Gradient.hpp"

Gradient::Gradient () {
  config_ = NULL;
  mesh_ = NULL;
  variableGradient_ = NULL;
}

Gradient::~Gradient () {
  if (variableGradient_ != NULL) {
    for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      for (short iVar=0; iVar<nVariable_; iVar++){
        delete [] variableGradient_[iElement][iVar];
      }
      delete [] variableGradient_[iElement];
    }
    delete [] variableGradient_;
  }
}

void Gradient::initialize (Config* config, Mesh* mesh, VariableType variableType) {
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
  
  variableGradient_ = new double** [ mesh_->getNumElement() ];
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    variableGradient_[iElement] = new double* [ nVariable_ ];
    for (short iVar=0; iVar<nVariable_; iVar++) {
      variableGradient_[iElement][iVar] = new double [ mesh_->getNumDim() ];
      for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
        variableGradient_[iElement][iVar][iDim] = 0;
      }
    }
  }
}

double Gradient::getGradient(
  const unsigned long &iElement, const short &iVariable, const short &iDim) {
  
  return variableGradient_[iElement][iVariable][iDim];
}
