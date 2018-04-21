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

#include "../include/SolutionSnapshot.hpp"

SolutionSnapshot::SolutionSnapshot() {
  variableContainer_ = NULL;
  conservedVariables_ = NULL;
  primitiveVariables_ = NULL;
  dimensionalPrimitiveVariables_ = NULL;
  config_ = NULL;
  gas_ = NULL;
  mesh_ = NULL;
}

SolutionSnapshot::~SolutionSnapshot() {
  
  if (conservedVariables_ != NULL) {
    for (short iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      delete [] conservedVariables_[iElement];
    }
  }
  delete [] conservedVariables_;
  
  if (primitiveVariables_ != NULL) {
    for (short iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      delete [] primitiveVariables_[iElement];
    }
  }
  delete [] primitiveVariables_;
  
  if (dimensionalPrimitiveVariables_ != NULL) {
    for (short iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      delete [] dimensionalPrimitiveVariables_[iElement];
    }
  }
  delete [] dimensionalPrimitiveVariables_;
  
  if (variableContainer_ != NULL) {
    delete [] variableContainer_;
  }
  
}

void SolutionSnapshot::initialize(Config *config, Mesh *mesh, Gas *gas) {
  
  config_ = config;
  gas_ = gas;
  mesh_ = mesh;
  
  conservedVariables_ = new double* [ mesh_->getNumElement() ];
  primitiveVariables_ = new double* [ mesh_->getNumElement() ];
  dimensionalPrimitiveVariables_ = new double* [ mesh_->getNumElement() ];
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    conservedVariables_[iElement] = new double [ config_->getNumConservedVariable() ];
    primitiveVariables_[iElement] = new double [ config_->getNumPrimitiveVariable() ];
    dimensionalPrimitiveVariables_[iElement] = new double [ config_->getNumPrimitiveVariable() ];
    for (short iVar=0; iVar<config_->getNumConservedVariable(); iVar++) {
      conservedVariables_[iElement][iVar] = 0;
    }
    for (short iVar=0; iVar<config_->getNumPrimitiveVariable(); iVar++) {
      primitiveVariables_[iElement][iVar] = 0;
      dimensionalPrimitiveVariables_[iElement][iVar] = 0;
    }
  }
  
  variableContainer_ = new double** [VariableType::COUNT];
  variableContainer_[VariableType::CONSERVED] = conservedVariables_;
  variableContainer_[VariableType::PRIMITIVE] = primitiveVariables_;
  
}

void SolutionSnapshot::setConservedVariable (
  const unsigned long &iElement, const short &iVariable, double value) {
  conservedVariables_[iElement][iVariable] = value;
}

void SolutionSnapshot::setConservedVariable (
  const unsigned long &iElement, double *conservedVariable) {
  for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
    conservedVariables_[iElement][iVariable] = conservedVariable[iVariable];
  }
}

void SolutionSnapshot::computeAllPrimitiveVariable () {
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    gas_->computeAllPrimitiveVariable(
      conservedVariables_[iElement], primitiveVariables_[iElement]);
  }
}

double* SolutionSnapshot::getConservedVariable(const unsigned long &iElement) {
  return conservedVariables_[iElement];
}

double SolutionSnapshot::getConservedVariable(
  const unsigned long &iElement, const short &iVariable) {
  return conservedVariables_[iElement][iVariable];
}

double SolutionSnapshot::getConservedVariable(
  const unsigned long &iElement, FlowVariable flowVariable) {
  return conservedVariables_[iElement][config_->mapConservedVariable(flowVariable)];
}

void SolutionSnapshot::addToConservedVariable (
  const unsigned long &iElement, const short &iVariable, double value) {
  conservedVariables_[iElement][iVariable] += value;
}

double* SolutionSnapshot::getPrimitiveVariable(const unsigned long &iElement) {
  return primitiveVariables_[iElement];
}

double SolutionSnapshot::getPrimitiveVariable(
  const unsigned long &iElement, const short &iVariable) {
  return primitiveVariables_[iElement][iVariable];
}

double SolutionSnapshot::getPrimitiveVariable(
  const unsigned long &iElement, FlowVariable flowVariable) {
  return primitiveVariables_[iElement][config_->mapPrimitiveVariable(flowVariable)];
}

double SolutionSnapshot::getVariable(
  const unsigned long &iElement, const short &iVariable, VariableType variableType) {
  return variableContainer_[variableType][iElement][iVariable];
}

void SolutionSnapshot::writePrimitiveToCsv(string fileName) {
  
  ofstream outputFile;
  
  outputFile.open(fileName);
  if (outputFile) {
    for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      
      gas_->convertDimensional(primitiveVariables_[iElement], dimensionalPrimitiveVariables_[iElement]);
      
      for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
        outputFile << mesh_->getElement(iElement)->getCentroid()->getSingleCoord(iDim) << ",";
      }
      for (short iVariable=0; iVariable<config_->getNumPrimitiveVariable(); iVariable++) {
        outputFile << dimensionalPrimitiveVariables_[iElement][iVariable] << ",";
      }
      outputFile << endl;
    }
    outputFile.close();
  }
  else {
    cout << endl << "WARNING: Failed to write solution to \"" << fileName << "\".";
  }
  
  
}

void SolutionSnapshot::writePrimitiveToVtk(string fileName) {
  
  unsigned long nElementVertexInfo;
  
  ofstream outputFile;
  
  Element *tmpElement;
  
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    gas_->convertDimensional(primitiveVariables_[iElement], dimensionalPrimitiveVariables_[iElement]);
  }
  
  outputFile.open(fileName);
  if (outputFile) {
    
    // Vertex information
    outputFile << "# vtk DataFile Version 2.0" << endl;
    outputFile << fileName << endl;
    outputFile << "ASCII" << endl;
    outputFile << "DATASET UNSTRUCTURED_GRID" << endl;
    outputFile << "POINTS " << mesh_->getNumVertex() << " float" << endl;
    for (unsigned long iVertex=0; iVertex<mesh_->getNumVertex(); iVertex++) {
      for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
        outputFile << mesh_->getVertex(iVertex)->getSingleCoord(iDim) << " ";
      }
      if (mesh_->getNumDim() == 2) {
        outputFile << 0.0 << " ";
      }
      outputFile << endl;
    }
    
    // Element connectivity
    nElementVertexInfo = 0;
    for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      nElementVertexInfo += mesh_->getElement(iElement)->getNumVertex() + 1;
    }
    outputFile << "CELLS " << mesh_->getNumElement() << " " << nElementVertexInfo << endl;
    for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      tmpElement = mesh_->getElement(iElement);
      outputFile << tmpElement->getNumVertex() << " ";
      for (int iVertex=0; iVertex<tmpElement->getNumVertex(); iVertex++) {
        outputFile << tmpElement->getSingleVertex(iVertex)->getId() << " ";
      }
      outputFile << endl;
    }
    
    // Element type
    outputFile << "CELL_TYPES " << mesh_->getNumElement() << endl;
    for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      outputFile << mesh_->getElement(iElement)->getVtkType() << endl;
    }
    
    // Element data
    outputFile << "CELL_DATA " << mesh_->getNumElement() << endl;
    outputFile << "SCALARS rho_nondim double" << endl;
    outputFile << "LOOKUP_TABLE default" << endl;
    for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      outputFile << 
        dimensionalPrimitiveVariables_[iElement][config_->mapPrimitiveVariable(FlowVariable::RHO)] << endl;
    }
    
    outputFile << "SCALARS P_nondim double" << endl;
    outputFile << "LOOKUP_TABLE default" << endl;
    for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      outputFile << 
        dimensionalPrimitiveVariables_[iElement][config_->mapPrimitiveVariable(FlowVariable::P)] << endl;
    }
    
    outputFile << "SCALARS T_nondim double" << endl;
    outputFile << "LOOKUP_TABLE default" << endl;
    for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      outputFile << 
      dimensionalPrimitiveVariables_[iElement][config_->mapPrimitiveVariable(FlowVariable::T)] << endl;
    }
    
    outputFile << "VECTORS M double" << endl;
    for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
        outputFile << 
          dimensionalPrimitiveVariables_[iElement][config_->mapPrimitiveVelocityComponent(iDim)]
          << " ";
      }
      if (mesh_->getNumDim() == 2) {
        outputFile << 0.0 << " ";
      }
      outputFile << endl;
    }
    
    outputFile.close();
  }
  else {
    cout << endl << "WARNING: Failed to write solution to \"" << fileName << "\".";
  }
  
}

void SolutionSnapshot::writeSlipWallToCsv(string fileName) {
  // Zeroth order interpolation of variables to surface
  
  ofstream outputFile;
  
  Face* face;
  unsigned long elementId;
  
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    gas_->convertDimensional(primitiveVariables_[iElement], dimensionalPrimitiveVariables_[iElement]);
  }
  
  outputFile.open(fileName);
  if (outputFile) {
    
    for (unsigned long iFace=0; iFace<mesh_->getNumFace(); iFace++) {
      
      face = mesh_->getFace(iFace);
      if (face->getBoundaryType() == BoundaryType::SLIP_WALL) {
        elementId = face->getOwnerElement()->getId();
        for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
          outputFile << face->getCenter()->getSingleCoord(iDim) << ",";
        }
        for (short iDim=0; iDim<mesh_->getNumDim(); iDim++) {
          outputFile << face->getDimNormal(iDim) << ",";
        }
        for (short iVariable=0; iVariable<config_->getNumPrimitiveVariable(); iVariable++) {
          outputFile << dimensionalPrimitiveVariables_[elementId][iVariable] << ",";
        }
        outputFile << endl;
      } // if face->getBoundary() == SLIP_WALL
      
    } // for face
    
    outputFile.close();
    
  }
  else {
    cout << endl << "WARNING: Failed to write slip wall variables to \"" << fileName << "\".";
  }
  
}