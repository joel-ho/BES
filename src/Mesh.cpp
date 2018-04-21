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

#include "../include/Mesh.hpp"
#include "../include/MeshGmsh.hpp"

Mesh::Mesh () {
  nDim_ = 0;
  nActiveBoundary_ = 0;
  nVertex_ = 0;
  nFace_ = 0;
  nElement_ = 0;
  activeBoundaryType_ = NULL;
  vertices_ = NULL;
  faces_ = NULL;
  elements_ = NULL;
}
Mesh::~Mesh () {
  
  if (activeBoundaryType_ != NULL) { 
    delete [] activeBoundaryType_;
    activeBoundaryType_ = NULL;
  }
  
  if (vertices_ != NULL) { 
    for (unsigned long iVertex=0; iVertex<nVertex_; iVertex++) {
      delete vertices_[iVertex];
    }
    delete [] vertices_;
    vertices_ = NULL;
  }
  
  if (faces_ != NULL) { 
    for (unsigned long iFace=0; iFace<nFace_; iFace++) {
      delete faces_[iFace];
    }
    delete [] faces_;
    faces_ = NULL;
  }
  
  if (elements_ != NULL) { 
    for (unsigned long iElement=0; iElement<nElement_; iElement++) {
      delete elements_[iElement];
    }
    delete [] elements_;
    elements_ = NULL;
  }
  
}

void Mesh::initialize (Config *config) {
  
  switch ( config->getMeshFormat() ) {
    
    case MeshFormat::GMSH:
      cout << "Initializing mesh from Gmsh file: \""<< config->getMeshFilePath() <<"\"... " << endl;
      readGmsh( config->getMeshFilePath() );
      break;
      
    default:
      cout << "ERROR: Mesh format not supported. Exiting..." << endl;
      exit(EXIT_FAILURE);
      break;
    
  }
  
  completeInitialization_();
  
  checkMesh_();
  
  cout << "  Number of dimensions: " << nDim_ << endl;
  cout << "  Number of vertices: " << nVertex_ << endl;
  cout << "  Number of faces: " << nFace_ << endl;
  cout << "  Number of elements: " << nElement_ << endl;
  
  cout << "Mesh initialized." << endl << endl;

}

void Mesh::readGmsh (string meshFilePath) {
  // Function sets:
  //    nDim_
  //    nVertex_
  //    nFace_
  //    nElement_
  //    vertices_
  //    faces_, computeCenter, computeNormal
  //    elements_
  //    Sets correct orientation of face normal (point from owner to neighbor)
  // Nodes to be entered into element in vtk ordering
  
  short numTags, tmpNumDim;
  int nBoundary=0, elementType;
  unsigned long i=0, n=-1;
  double *elementOutwardNormal, *faceNormal;
  string fileLine, faceIdentifier;
  string physicalHeader, nodesHeader, elementsHeader;
  string delimiter=" ";
  
  string splitLine[LINE_BUFFER];
  double splitLineNum[LINE_BUFFER];
  for (int iEntry=0; iEntry<LINE_BUFFER; iEntry++) {
    splitLine[iEntry] = "";
    splitLineNum[iEntry] = 0;
  }
  
  ifstream mshFile;
  
  Face *face;
  Element *element;
  
  vector <unsigned long> faceVertexIndices;
  set <BoundaryType> activeBoundaryType;
  vector <Face*> mshFileFaces;
  vector <Element*> mshFileElements;
  
  enum Section { NIL, PHYSICAL_NAMES, NODES, ELEMENTS };
  Section currentSection = NIL;
  
  // Maps gmsh boundary ID to code's boundary ID through map defined in Options.cpp
  map <int, BoundaryType> gmshBoundaryTag; 
  
  // Maps face identifier to face pointer
  map <string, Face*> gmshFaceTag; 
  
  // Element face node ordering for outward normal
  short gmshFourNodeOrdering[] = {
    0, 2, 1, 
    1, 2, 3, 
    0, 3, 2, 
    0, 1, 3};
  
  physicalHeader = "$PhysicalNames";
  nodesHeader = "$Nodes";
  elementsHeader = "$Elements";
  
  mshFile.open(meshFilePath);
  while ( getline(mshFile, fileLine) ) {
    
    if ( fileLine.compare(physicalHeader) == 0 ) {
      currentSection = PHYSICAL_NAMES;
      n = 0;
      continue;
    }
    else if ( fileLine.compare(nodesHeader) == 0 ) {
      currentSection = NODES;
      n = 0;
      continue;
    }
    else if ( fileLine.compare(elementsHeader) == 0 ) {
      currentSection = ELEMENTS;
      n = 0;
      continue;
    }
    
    if ( currentSection == PHYSICAL_NAMES ) {
      
      if ( n == 0 ) {
        n = stod(fileLine);
        nBoundary = n;
      }
      else if ( i < n ) {
        
        ReadFile::splitString(fileLine, delimiter, splitLine);
        
        tmpNumDim = (short) stod(splitLine[0]);
        nDim_ = ( tmpNumDim > nDim_ )? tmpNumDim : nDim_;
        
        gmshBoundaryTag[ (int) stoi(splitLine[1]) ] = BoundaryTypeStringMap[
          splitLine[2].substr(1, splitLine[2].size()-2)];
        i++;
        
      }
      else {
        currentSection = NIL;
        n = -1;
        i = 0;
      }
    }
    
    if ( currentSection == NODES ) {
      
      if ( n == 0 ) {
        n = stod(fileLine);
        nVertex_ = n;
        vertices_ = new Vertex* [nVertex_];
      }
      else if ( i < n ) {
        
        ReadFile::splitStringToNum(fileLine, delimiter, splitLineNum);
        
        vertices_[i] = new Vertex();
        vertices_[i]->initialize(nDim_);
        for ( short iDim=0; iDim<nDim_; iDim++ ) {
          vertices_[i]->setSingleCoord(splitLineNum[iDim+1], iDim);
        }
        vertices_[i]->setId(i);
        i++;
        
      }
      else {
        currentSection = NIL;
        n = -1;
        i = 0;
      }
    }
    
    if ( currentSection == ELEMENTS ) {
      
      if ( n == 0 ) {
        n = stod(fileLine);
        mshFileElements.reserve(n);
        mshFileFaces.reserve(n);
        elementOutwardNormal = new double [nDim_];
      }
      else if ( i < n ) {
        ReadFile::splitStringToNum(fileLine, delimiter, splitLineNum);
        elementType = splitLineNum[1];
        numTags = splitLineNum[2];
        
        switch (elementType) {
          
          case 1: // Line element - line is face in 2D case
            
            if (nDim_ > 2) {
              cout << "ERROR: 1D boundary encountered in 3D problem. Exiting..." << endl;
              exit(EXIT_FAILURE);
            }
            
            // Extract vertices number
            faceVertexIndices.clear();
            for (short m=0; m<2; m++) {
              faceVertexIndices.push_back(
                (unsigned long) splitLineNum[3+numTags+m]-1 );
            }
            
            // Create new Face object
            face = new FaceGmshOne();
            mshFileFaces.push_back(face);
            mshFileFaces.back()->initialize();
            for (short m=0; m<2; m++) {
              mshFileFaces.back()->setSingleVertex(
                vertices_[ faceVertexIndices[m] ], m);
            }
              
            mshFileFaces.back()->setBoundaryType(
              gmshBoundaryTag[ (int) splitLineNum[1+numTags] ] );
            
            mshFileFaces.back()->computeCenter();
            mshFileFaces.back()->computeNormal();
            
            // Create set of all active boundary types in mesh
            activeBoundaryType.insert(
              gmshBoundaryTag[ (int) splitLineNum[1+numTags] ] );
            
            // Add to log of all faces to create face based data storage system
            ReadFile::encodeFaceIdentifier(faceIdentifier, faceVertexIndices);
            gmshFaceTag[faceIdentifier] = mshFileFaces.back();
            
            break;
            
          case 2: // Triangle element
          
            if (nDim_ > 2) { // 3D case - triangle is face
            
              // Extract vertices number
              faceVertexIndices.clear();
              for (short m=0; m<3; m++) {
                faceVertexIndices.push_back(
                  (unsigned long) splitLineNum[3+numTags+m]-1 );
              }
              
              // Create new Face object
              face = new FaceGmshTwo();
              mshFileFaces.push_back(face);
              mshFileFaces.back()->initialize();
              for (short m=0; m<3; m++) {
                mshFileFaces.back()->setSingleVertex(
                  vertices_[ faceVertexIndices[m] ], m);
              }
              mshFileFaces.back()->setBoundaryType(
                gmshBoundaryTag[ (int) splitLineNum[1+numTags] ] );
              
              mshFileFaces.back()->computeCenter();
              mshFileFaces.back()->computeNormal();
              
              // Create set of all active boundary types in mesh
              activeBoundaryType.insert(
                gmshBoundaryTag[ (int) splitLineNum[1+numTags] ] );
              
              // Add to log of all faces to create face based data storage system
              ReadFile::encodeFaceIdentifier(faceIdentifier, faceVertexIndices);
              gmshFaceTag[faceIdentifier] = mshFileFaces.back();
              
            }
            
            else { // 2D case - triangle is element
              
              // Add element to mesh
              element = new ElementGmshTwo();
              mshFileElements.push_back(element);
              
              // Add vertices to element
              // Vertices to be added in vtk element order
              mshFileElements.back()->initialize();
              for (short m=0; m<3; m++) {
                mshFileElements.back()->setSingleVertex(
                  vertices_[ (unsigned long)splitLineNum[3+numTags+m] - 1 ], m); // Subtract 1 from vertex number from msh file due to different starting index (1 for msh file, and 0 for array here)
              }
               
              // Add element faces to mesh faces.
              // Loop across all faces in element.
              for (short m=0; m<3; m++) {
                
                // Get face identifier
                faceVertexIndices.clear();
                faceVertexIndices.push_back(
                  (unsigned long)splitLineNum[3+numTags+m] - 1);
                if (m<2) {
                  faceVertexIndices.push_back(
                    (unsigned long)splitLineNum[3+numTags+m+1] - 1);
                }
                else {
                  faceVertexIndices.push_back(
                    (unsigned long)splitLineNum[3+numTags] - 1);
                }
                ReadFile::encodeFaceIdentifier(faceIdentifier, faceVertexIndices);
                
                // If face already exists
                if (gmshFaceTag.count(faceIdentifier) == 1) { 
                  
                  // Calculate element outward normal on particular face
                  if (m<2) {
                    MeshMath::computeTwoDimOutwardNormal(
                      elementOutwardNormal,
                      mshFileElements.back()->getSingleVertex(m),
                      mshFileElements.back()->getSingleVertex(m+1));
                  }
                  else {
                    MeshMath::computeTwoDimOutwardNormal(
                      elementOutwardNormal,
                      mshFileElements.back()->getSingleVertex(m),
                      mshFileElements.back()->getSingleVertex(0));
                  }
                  
                  // Get face normal
                  faceNormal = gmshFaceTag[faceIdentifier]->getDimNormal();
                  
                  // If face belongs to external boundary 
                  if (gmshFaceTag[faceIdentifier]->queryBoundaryExternal() == 1) { 
                    
                    // Check if element outward normal and existing face normal
                    // are aligned. Flip existing face normal otherwise.
                    if (MeshMath::dotProduct(faceNormal, elementOutwardNormal, nDim_) < 0) {
                      gmshFaceTag[faceIdentifier]->flipNormal();
                    }
                    
                    // Store pointer to element and set neighbor element to 
                    // NULL for external boundary 
                    gmshFaceTag[faceIdentifier]->setOwnerElement(mshFileElements.back());
                    gmshFaceTag[faceIdentifier]->setNeighborElement(NULL);
                    
                  }
                  
                  // If face belongs to internal boundary 
                  else if (gmshFaceTag[faceIdentifier]->queryBoundaryExternal() == -1) {
                    gmshFaceTag[faceIdentifier]->setNeighborElement(mshFileElements.back());
                  }
                  
                  // Boundary entered by user not supported
                  else {
                    cout << "ERROR: Unknown boundary type entered. Exiting..." << endl;
                    exit(EXIT_FAILURE);
                  }
                
                }
                
                // Face does not exist
                else {
                  
                  // Create new Face object
                  mshFileFaces.push_back(new FaceGmshOne);
                  mshFileFaces.back()->initialize();
                  
                  // Enter vertex pointers (take note of order)
                  if (m<2) {
                    mshFileFaces.back()->setSingleVertex(
                      mshFileElements.back()->getSingleVertex(m+1), 1);
                  }
                  else {
                    mshFileFaces.back()->setSingleVertex(
                      mshFileElements.back()->getSingleVertex(0), 1);
                  }
                  
                  mshFileFaces.back()->setSingleVertex(
                    mshFileElements.back()->getSingleVertex(m), 0);
                    
                  mshFileFaces.back()->setBoundaryType(
                    BoundaryType::INTERNAL );
                  mshFileFaces.back()->computeCenter();
                  mshFileFaces.back()->computeNormal(); // With order that vertices are entered, normal will be pointing outward
                  
                  // Add to log of all faces to create face based data storage system
                  ReadFile::encodeFaceIdentifier(faceIdentifier, faceVertexIndices);
                  gmshFaceTag[faceIdentifier] = mshFileFaces.back();
                  
                  // Add element pointer to face
                  gmshFaceTag[faceIdentifier]->setOwnerElement(mshFileElements.back());
                  
                }
              
                // Save face to element
                mshFileElements.back()->setSingleFace(gmshFaceTag[faceIdentifier], m);
              
              }
            }
            break;
            
          case 4: // Tetrahedron as element
            
            // Add element to mesh
            element = new ElementGmshFour();
            mshFileElements.push_back(element);
            
            // Add vertices to element
            // Vertices to be added in vtk element order
            mshFileElements.back()->initialize();
            for (short m=0; m<4; m++) {
              mshFileElements.back()->setSingleVertex(
                vertices_[ (unsigned long)splitLineNum[3+numTags+m] - 1 ], m); // Subtract 1 from vertex number from msh file due to different starting index (1 for msh file, and 0 for array here)
            }
            
            // Add element faces to mesh faces.
            // Loop across all faces in element.
            for (short m=0; m<4; m++) {
              
              // Get face identifier
              faceVertexIndices.clear();
              for (short iVertex=0; iVertex<3; iVertex++) {
                faceVertexIndices.push_back(
                (unsigned long)splitLineNum[ 3+numTags+gmshFourNodeOrdering[m*3+iVertex] ] - 1);
              }
              ReadFile::encodeFaceIdentifier(faceIdentifier, faceVertexIndices);
              
              // If face already exists
              if (gmshFaceTag.count(faceIdentifier) == 1) { 
                
                // Calculate element outward normal on particular face
                MeshMath::crossProduct(
                  elementOutwardNormal,
                  mshFileElements.back()->getSingleVertex(gmshFourNodeOrdering[m*3+0]),
                  mshFileElements.back()->getSingleVertex(gmshFourNodeOrdering[m*3+1]),
                  mshFileElements.back()->getSingleVertex(gmshFourNodeOrdering[m*3+2]));
                
                // Get face normal
                faceNormal = gmshFaceTag[faceIdentifier]->getDimNormal();
                
                // If face belongs to external boundary 
                if (gmshFaceTag[faceIdentifier]->queryBoundaryExternal() == 1) { 
                  
                  // Check if element outward normal and existing face normal
                  // are aligned. Flip existing face normal otherwise.
                  if (MeshMath::dotProduct(faceNormal, elementOutwardNormal, nDim_) < 0) {
                    gmshFaceTag[faceIdentifier]->flipNormal();
                  }
                  
                  // Store pointer to element and set neighbor element to 
                  // NULL for external boundary 
                  gmshFaceTag[faceIdentifier]->setOwnerElement(mshFileElements.back());
                  gmshFaceTag[faceIdentifier]->setNeighborElement(NULL);
                  
                }
                
                // If face belongs to internal boundary 
                else if (gmshFaceTag[faceIdentifier]->queryBoundaryExternal() == -1) {
                  gmshFaceTag[faceIdentifier]->setNeighborElement(mshFileElements.back());
                }
                
                // Boundary entered by user not supported
                else {
                  cout << "ERROR: Unknown boundary type entered. Exiting..." << endl;
                  exit(EXIT_FAILURE);
                }
              
              } // end if face already exists
              
              // Face does not exist
              else {
                
                // Create new Face object
                mshFileFaces.push_back(new FaceGmshTwo);
                mshFileFaces.back()->initialize();
                
                // Enter vertex pointers (take note of order)
                for (short iVertex=0; iVertex<3; iVertex++) {
                  mshFileFaces.back()->setSingleVertex(
                    mshFileElements.back()->getSingleVertex( gmshFourNodeOrdering[m*3+iVertex] ), iVertex);
                }
                
                mshFileFaces.back()->setBoundaryType(
                  BoundaryType::INTERNAL );
                mshFileFaces.back()->computeCenter();
                mshFileFaces.back()->computeNormal(); // With order that vertices are entered, normal will be pointing outward
                
                // Add to log of all faces to create face based data storage system
                ReadFile::encodeFaceIdentifier(faceIdentifier, faceVertexIndices);
                gmshFaceTag[faceIdentifier] = mshFileFaces.back();
                
                // Add element pointer to face
                gmshFaceTag[faceIdentifier]->setOwnerElement(mshFileElements.back());
                
              }
            
              // Save face to element
              mshFileElements.back()->setSingleFace(gmshFaceTag[faceIdentifier], m);
            
            } // end for all faces in element
            break;
        }
        i++;
      }
      else {
        currentSection = NIL;
        n = -1;
        i = 0;
      }
    
    } // end if ( currentSection == ELEMENTS )
    
  } // end while ( getline(mshFile, fileLine) )

  
  // // Gmsh msh file had element vertices entered in clockwise order,
  // // which is different from documentation and convention used above.
  // // Flip normals to correct for discrepency in node ordering.
  // for (unsigned long iFace=0; iFace<mshFileFaces.size(); iFace++) {
    // mshFileFaces[iFace]->flipNormal();
  // }
  
  saveFaceElement_(mshFileFaces, mshFileElements);
  saveActiveBoundary_(activeBoundaryType);
  
}

void Mesh::saveFaceElement_ (vector <Face*> &faces, vector <Element*> &elements) {
  nFace_ = faces.size();
  nElement_ = elements.size();
  faces_ = new Face* [nFace_];
  elements_ = new Element* [nElement_];
  copy(faces.begin(), faces.end(), faces_);
  copy(elements.begin(), elements.end(), elements_);
}

void Mesh::saveActiveBoundary_ (set<BoundaryType> &activeBoundaryType) {
  
  nActiveBoundary_ = activeBoundaryType.size();
  activeBoundaryType_ = new BoundaryType [nActiveBoundary_];
  for (
    auto boundaryIterator=activeBoundaryType.begin(); 
    boundaryIterator!=activeBoundaryType.end(); 
    boundaryIterator++) {
      
    auto iBoundary = distance(activeBoundaryType.begin(), boundaryIterator);
    activeBoundaryType_[iBoundary] = *boundaryIterator;
    
  }
  
}

void Mesh::completeInitialization_ () {
  
  unsigned long elementId;
  Face* elementFace;
  
  for (unsigned long iFace=0; iFace<nFace_; iFace++) {
    faces_[iFace]->setId(iFace);
  }
  
  for (unsigned long iElement=0; iElement<nElement_; iElement++) {
    elements_[iElement]->setId(iElement);
    elements_[iElement]->computeCentroid();
    elements_[iElement]->computeVolume();
    elements_[iElement]->computeAllVectorToFace();
  }
  
  for (unsigned long iElement=0; iElement<nElement_; iElement++) {
    
    // Set pointers to neighboring elements
    elementId = elements_[iElement]->getId();
    for (short iFace=0; iFace<elements_[iElement]->getNumFace(); iFace++) {
      elementFace = elements_[iElement]->getSingleFace(iFace);
      
      // Internal face
      if (elementFace->queryBoundaryExternal() == -1) {
        
        // If current element is not owner of face, 
        // set owner of face as adjacent element to face
        if (elementId != elementFace->getOwnerElement()->getId()) {
        elements_[iElement]->setFaceAdjacentElement(elementFace->getOwnerElement(), iFace);
        }
        
        // If current element is owner of face, 
        // set neighbor of face as adjacent element to face
        else {
          elements_[iElement]->setFaceAdjacentElement(elementFace->getNeighborElement(), iFace);
        }
      }
      
      // External face
      else if (elementFace->queryBoundaryExternal() == 1) {
        elements_[iElement]->setFaceAdjacentElement(NULL, iFace);
      }
      
    }
    
    // Calculate inverse distance weights to adjoining face
    elements_[iElement]->computeAllFaceInverseDistanceWeightage();
    
  }
  
}

void Mesh::checkMesh_() {
  
  if (nVertex_ == 0) {
    cout << "ERROR: No vertices in mesh file. Exiting..." << endl;
    exit(EXIT_FAILURE);
  }
  
  if (nFace_ == 0) {
    cout << "ERROR: No faces in mesh file. Exiting..." << endl;
    exit(EXIT_FAILURE);
  }
  
  if (nElement_ == 0) {
    cout << "ERROR: No elements in mesh file. Exiting..." << endl;
    exit(EXIT_FAILURE);
  }
  
}

short Mesh::getNumDim() {
  return nDim_;
}

short Mesh::getNumActiveBoundary () {
  return nActiveBoundary_;
}

BoundaryType Mesh::getActiveBoundaryType (const short &iBoundary) {
  return activeBoundaryType_[iBoundary];
}

Vertex* Mesh::getVertex(const unsigned long &iVertex) {
  return vertices_[iVertex];
}

Face* Mesh::getFace(const unsigned long &iFace) {
  return faces_[iFace];
}

Element* Mesh::getElement(const unsigned long &iElement) {
  return elements_[iElement];
}

void Mesh::getElement(Element *resultantElement, const unsigned long &iElement) {
  resultantElement = elements_[iElement];
}

unsigned long Mesh::getNumVertex() {
  return nVertex_;
}

unsigned long Mesh::getNumFace() {
  return nFace_;
}

unsigned long Mesh::getNumElement() {
  return nElement_;
}

void Mesh::printCoords() {
  for (unsigned long i=0; i<nVertex_; i++) {
    vertices_[i]->printCoords(nDim_);
  }
}
