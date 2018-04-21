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

#include "../include/SolverImplicit.hpp"

//--------------------------//
// parfor private variables //
//--------------------------//
FaceLoopVariables::FaceLoopVariables() {}

FaceLoopVariables::~FaceLoopVariables() {
  for 
  (
    map <BoundaryType, Flux*>::iterator it = fluxCalculator.begin();
    it != fluxCalculator.end();
    ++it
  ) 
  {
    delete it->second;
    delete fluxJacobianCalculator[it->first];
  }
}

void FaceLoopVariables::initialize
(
  Config* configIn, 
  Gas* gasIn, 
  Mesh* meshIn, 
  unsigned long iTimeStepIn,
  int printOptions
)
{
  
  BoundaryType activeBoundaryType;
  
  // Save properties
  config = configIn;
  gas = gasIn;
  mesh = meshIn;
  iTimeStep = iTimeStepIn;
  
  // Select required boundary flux and flux jacobian calculators
  for (short iBoundary=0; iBoundary<mesh->getNumActiveBoundary(); iBoundary++) {
    
    activeBoundaryType = mesh->getActiveBoundaryType(iBoundary);
    
    switch (activeBoundaryType) {  
      
      case BoundaryType::FREESTREAM:
        
        if (printOptions > 0) {
          cout << "  Initializing freestream boundary... ";
        }
        
        fluxCalculator[activeBoundaryType] = new FluxFreestream;
        rightVariableMap[activeBoundaryType] = 
          config->getFreestreamPrimitiveVariables(iTimeStep);
          
        if (config->getNumDim() == 2) {
          fluxJacobianCalculator[activeBoundaryType] = new FluxJacobianFreestreamTwoDim;
        }
        else if (config->getNumDim() == 3) {
          fluxJacobianCalculator[activeBoundaryType] = new FluxJacobianFreestreamThreeDim;
        }
        rightJacobianMap[activeBoundaryType] = NULL;
        
        if (printOptions > 0) {
          cout << "Done." << endl;
        }
        
        break;
      
        
      case BoundaryType::SLIP_WALL:
      
        if (printOptions > 0) {
          cout << "  Initializing slip wall boundary... ";
        }
      
        fluxCalculator[activeBoundaryType] = new FluxSlipWall;
        rightVariableMap[activeBoundaryType] = NULL;
        
        if (config->getNumDim() == 2) {
          fluxJacobianCalculator[activeBoundaryType] = new FluxJacobianSlipWallTwoDim;
        }
        else if (config->getNumDim() == 3) {
          fluxJacobianCalculator[activeBoundaryType] = new FluxJacobianSlipWallThreeDim;
        }
        rightJacobianMap[activeBoundaryType] = NULL;
        
        if (printOptions > 0) {
          cout << "Done." << endl;
        }
        
        break;
    
    }
  }
  
  // Select internal flux and flux jacobian calculator
  switch (config->getFluxDiscretizationScheme()) {
  
    case FluxDiscretizationScheme::ROE:
    
      if (printOptions > 0) {
        cout << "  Initializing Roe flux difference splitting... ";
      }
      
      fluxCalculator[BoundaryType::INTERNAL] = new FluxRoe;
      rightVariableMap[BoundaryType::INTERNAL] = rightVariableBuffer;
      
      if (config->getNumDim() == 2) {
        fluxJacobianCalculator[BoundaryType::INTERNAL] = new FluxJacobianRoeTwoDim;
      }
      else if (config->getNumDim() == 3) {
        fluxJacobianCalculator[BoundaryType::INTERNAL] = new FluxJacobianRoeThreeDim;
      }
      rightJacobianMap[BoundaryType::INTERNAL] = rightJacobianBuffer;
            
      if (printOptions > 0) {
        cout << "Done." << endl;
      }
      
      break;
      
  }
  
  // Initialize flux and flux jacobian calculators
  for 
  (
    map <BoundaryType, Flux*>::iterator it = fluxCalculator.begin();
    it != fluxCalculator.end();
    ++it
  ) 
  {
    it->second->initialize(config, gas);
    fluxJacobianCalculator[it->first]->initialize(config);
  }
  
}

FaceLoopVariables::FaceLoopVariables
(
  const FaceLoopVariables& solverImplicitVariables
) 
{
  initialize(
    solverImplicitVariables.config, 
    solverImplicitVariables.gas,
    solverImplicitVariables.mesh,
    solverImplicitVariables.iTimeStep,
    -1
  );
}

void FaceLoopVariables::printJacobian(short nVar, double* jacobian) {
  for (short i=0; i<nVar; i++) {
    for (short j=0; j<nVar; j++) {
      cout << jacobian[ i*nVar + j ] << ", ";
    }
    cout << endl;
  }
}

//-----------------//
// Implicit solver //
//-----------------//
void SolverImplicit::initialize (Config *config, Gas *gas, Mesh *mesh) {
  
  cout << "Initializing implicit solver... " << endl;
  
  BoundaryType activeBoundaryType;
  
  /*Solver->*/initializeBase_(config, gas, mesh);
  
  coefficientMatrix_ = new CoefficientMatrix;
  coefficientMatrix_->initialize(config_, mesh_);
  
  faceLoopVariables_.initialize(config, gas, mesh, 0, 1); // configIn, gasIn, meshIn, iTimeStepIn, printOptions
  
  switch (config_->getGradientScheme()) {
  
    case GradientScheme::GREEN_GAUSS:
      cout << "  Initializing Green-Gauss gradient... ";
      gradient_ = new GradientGreenGauss;
      break;
  
  }
  gradient_->initialize(config_, mesh_, VariableType::CONSERVED);
  cout << "Done." << endl;
  
  switch (config_->getLimiterScheme()) {
  
    case LimiterScheme::VENKATAKRISHNAN:
      cout << "  Initializing Venkatakrishnan limiter... ";
      limiter_ = new LimiterVenkatakrishnan;
      break;
  
  }
  limiter_->initialize(config_, mesh_, VariableType::CONSERVED);
  cout << "Done." << endl;

  switch (config_->getLocalTimeStepOption()) {
  
    case LocalTimeStep::YES:
      cout << "  Initializing local time stepper ";
      timeStepCalculator_ = new TimeStepCalculatorLocal;
      break;
      
    case LocalTimeStep::NO:
      cout << "  Initializing global time stepper ";
      timeStepCalculator_ = new TimeStepCalculatorGlobal;
      break;
  
  }
  timeStepCalculator_->initialize(config_, gas_, mesh_);
  
  if (config_->getAdaptiveCourantNumOption() == AdaptiveCourantNum::YES) {
    cout << "with adaptive Courant number..." << endl;
    cout << "    Min Courant number: " << config_->getCourantNum() << "..." << endl;
    cout << "    Max Courant number: " << config_->getAdaptiveCourantNumMax() << "..." << endl;
    cout << "    Increment exponent: " << config_->getAdaptiveCourantNumExponent()[0] << "..." << endl;
    cout << "    Decrement exponent: " << config_->getAdaptiveCourantNumExponent()[1] << "... ";
  }
  else {
    cout << "with Courant number: " << config_->getCourantNum() << "... ";
  }
  cout << "Done." << endl;
  
  explicitFlux_.resize(
    mesh->getNumElement()*config->getNumConservedVariable(), 0);
  
  deltaConservedVariable_.resize(
    mesh->getNumElement()*config->getNumConservedVariable(), 0);
    
  maxResidual_ = new double [config->getNumConservedVariable()];
  for (short iVar=0; iVar<config_->getNumConservedVariable(); iVar++) {
    maxResidual_[iVar] = 0;
  }
  maxResidualOld_ = new double [config->getNumConservedVariable()];
  
  actualCourantNum_ = config_->getCourantNum();
  
  cout << "Solver initialized." << endl << endl;
  
}

SolverImplicit::~SolverImplicit () {
  delete coefficientMatrix_;
  delete gradient_;
  delete limiter_;
  delete timeStepCalculator_;
  delete [] maxResidual_;
  delete [] maxResidualOld_;
}

void SolverImplicit::setUniversalInitialCondition(double *primitiveVariable) {
  
  cout << "Initializing flow variables in entire domain... " << endl;
  
  double conservedVariable[VAR_BUFFER];
  
  gas_->computeAllConservedVariable(primitiveVariable, conservedVariable);
  
  #pragma omp parallel for default(shared)
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    solution_->setConservedVariable(0, 0, iElement, conservedVariable);
  }
  
  cout << "Flow variables initialized." << endl << endl;
  
}

unsigned long SolverImplicit::computeVectorIdx_(unsigned long iElement, short iVar) {
  return iElement*config_->getNumConservedVariable() + iVar;
}

void SolverImplicit::timeStep_ 
(
  const unsigned long &iTimeStepStart,
  FaceLoopVariables loopVar,
  string residualFileName
) 
{
  
  int solverIter;
  double solverError;
  
  ostringstream residualStr;
  SolutionSnapshot *solutionSnapshot;
  
  resetSystemOfEquation_();
  
  solutionSnapshot = solution_->getCurrentSnapshot(iTimeStepStart, 0);
  
  timeStepCalculator_->computeAllTimeStep(
    solutionSnapshot,
    actualCourantNum_
  );
    
  gradient_->computeGradient(solutionSnapshot);
  
  limiter_->computeLimiterTerm(solutionSnapshot, gradient_);
  
  computeFaceFluxes_(iTimeStepStart, solutionSnapshot, loopVar);
  
  coefficientMatrix_->solve(explicitFlux_, deltaConservedVariable_, solverIter, solverError);
  
  updateSolution_(solutionSnapshot);
  
  adaptCourantNum_();
  
  // Print max residual and solver errors
  for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
    cout << scientific << maxResidual_[iVariable] << " | ";
    residualStr << scientific << maxResidual_[iVariable] << ",";
  }
  
  cout << solverIter << " | ";
  residualStr << solverIter << ",";
  
  cout << scientific << solverError << " | ";
  residualStr << scientific << solverError << ",";
  
  cout << endl;
  residualStr << endl;
  
  /*Solver->*/appendToResidual_(residualFileName, residualStr.str());
  
}

void SolverImplicit::resetSystemOfEquation_() {
  fill(explicitFlux_.begin(), explicitFlux_.end(), 0);
  fill(deltaConservedVariable_.begin(), deltaConservedVariable_.end(), 0);
  coefficientMatrix_->changeAllVal(0);
  coefficientMatrix_->addToDiagonal(1); // Identity matrix
  for (short iVar=0; iVar<config_->getNumConservedVariable(); iVar++) {
    maxResidualOld_[iVar] = maxResidual_[iVar];
    maxResidual_[iVar] = 0;
  }
}

void SolverImplicit::computeFaceFluxes_
(
  const unsigned long & iTimeStepStart,
  SolutionSnapshot * solutionSnapshot,
  FaceLoopVariables loopVar
) 
{
  
  loopVar.iTimeStep = iTimeStepStart;
  
  #pragma omp parallel for default(shared) firstprivate(loopVar)
  for (unsigned long iFace=0; iFace<mesh_->getNumFace(); iFace++) {
    
    // deltaT/volume*faceArea
    double multiplierLeft, multiplierRight;
    
    // Get properties at face
    loopVar.face = mesh_->getFace(iFace);
    loopVar.boundaryType = loopVar.face->getBoundaryType();
    loopVar.faceArea = loopVar.face->getArea();
    loopVar.leftElement = loopVar.face->getOwnerElement();
    loopVar.leftElementId = loopVar.leftElement->getId();
    loopVar.rightElement = loopVar.face->getNeighborElement();
    if (loopVar.rightElement != NULL) {
      loopVar.rightElementId = loopVar.rightElement->getId();
    }
    
    // Extrapolate variables to face, output to leftVariable and rightVariable
    /*Solution->*/computeFaceConservedVariable_(
      solutionSnapshot, 
      gradient_, 
      limiter_, 
      iFace, 
      loopVar.leftVariableBuffer, 
      loopVar.rightVariableBuffer);
    
    // Compute explicit flux across face, output to explicit flux buffer
    loopVar.fluxCalculator[loopVar.boundaryType]->computeFlux(
      loopVar.leftVariableBuffer, 
      loopVar.rightVariableMap[loopVar.boundaryType], 
      loopVar.face->getNonDimNormal(), 
      loopVar.fluxBuffer
    );
    
    // Get first order variables on face
    /*Solver->*/getFirstOrderFaceVariable_(
      solutionSnapshot, 
      iFace, 
      loopVar.leftVariableBuffer, 
      loopVar.rightVariableBuffer
    );
    
    // Compute first order jacobian
    loopVar.fluxJacobianCalculator[loopVar.boundaryType]->computeJacobian(
      loopVar.leftVariableBuffer, 
      loopVar.rightVariableMap[loopVar.boundaryType],
      loopVar.face->getNonDimNormal(),
      loopVar.leftJacobianBuffer,
      loopVar.rightJacobianMap[loopVar.boundaryType]
    );
    
    // Update RHS and coefficient matrix
    multiplierLeft = 
      timeStepCalculator_->getTimeStep(loopVar.leftElementId)/
      loopVar.leftElement->getVolume()*
      loopVar.faceArea;
      
    for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
      
      // Check for divergence
      if (loopVar.fluxBuffer[iVariable] != loopVar.fluxBuffer[iVariable]) {
        #pragma omp critical 
        {
        cout << endl << "ERROR: Nan values encountered in flux across face at coordinates [";
        for (short iDimError=0; iDimError<mesh_->getNumDim()-1; iDimError++) {
          cout << scientific << loopVar.face->getCenter()->getSingleCoord(iDimError) << ", ";
        }
        cout << scientific << loopVar.face->getCenter()->getSingleCoord(mesh_->getNumDim()-1);
        cout << "]. Exiting..." << endl;
        exit(EXIT_FAILURE);
        }
      }
      
      // Update element explicit flux (iOwner)
      #pragma omp atomic
      explicitFlux_[ computeVectorIdx_(loopVar.leftElementId, iVariable) ] -= 
        multiplierLeft*loopVar.fluxBuffer[iVariable];
        
      // Update coefficient matrix (iOwner, iOwner)
      // #pragma omp atomic directive in MatrixSparse class
      coefficientMatrix_->addToBlock(
        loopVar.leftElementId,
        loopVar.leftElementId,
        multiplierLeft,
        loopVar.leftJacobianBuffer
      );
      
    } // for iVariable
    
    // Update RHS and coefficient matrix for internal face
    if (loopVar.rightElement != NULL) {
      
      multiplierRight = 
        timeStepCalculator_->getTimeStep(loopVar.rightElementId)/
        loopVar.rightElement->getVolume()*
        loopVar.faceArea;
      
      for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
        
        // Update element explicit flux (iNeighbor)
        #pragma omp atomic
        explicitFlux_[ computeVectorIdx_(loopVar.rightElementId, iVariable) ] += 
          timeStepCalculator_->getTimeStep(loopVar.rightElementId)/
          loopVar.rightElement->getVolume()*
          loopVar.faceArea*
          loopVar.fluxBuffer[iVariable];
          
        // Update coefficient matrix (iOwner, iNeightbor)
        // #pragma omp atomic directive in MatrixSparse class
        coefficientMatrix_->addToBlock(
          loopVar.leftElementId,
          loopVar.rightElementId,
          multiplierLeft,
          loopVar.rightJacobianBuffer
        );
        
        // Update coefficient matrix (iNeightbor, iOwner)
        // #pragma omp atomic directive in MatrixSparse class
        coefficientMatrix_->addToBlock(
          loopVar.rightElementId,
          loopVar.leftElementId,
          -multiplierRight,
          loopVar.leftJacobianBuffer
        );
          
        // Update coefficient matrix (iNeightbor, iNeightbor)
        // #pragma omp atomic directive in MatrixSparse class
        coefficientMatrix_->addToBlock(
          loopVar.rightElementId,
          loopVar.rightElementId,
          -multiplierRight,
          loopVar.rightJacobianBuffer
        );
          
      } // for iVariable
      
    } // if loopVar.rightElement != NULL
    
  } // #pragma omp parallel for iFace
}

void SolverImplicit::updateSolution_(SolutionSnapshot * solutionSnapshot) {
  
  #pragma omp parallel for default(shared)
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {      
      
      solutionSnapshot->addToConservedVariable(
        iElement,
        iVariable,
        deltaConservedVariable_[ computeVectorIdx_(iElement, iVariable) ] );
        
      #pragma omp critical
      maxResidual_[iVariable] = 
        ( abs(explicitFlux_[ computeVectorIdx_(iElement, iVariable) ]) > maxResidual_[iVariable] ) ?
        abs(explicitFlux_[ computeVectorIdx_(iElement, iVariable) ]) : maxResidual_[iVariable];
        
    }
  }
  
}

void SolverImplicit::adaptCourantNum_() {
  
  short rhoIdx, exponentIdx;
  
  if (config_->getAdaptiveCourantNumOption() == AdaptiveCourantNum::YES) {
    
    rhoIdx = config_->mapConservedVariable(FlowVariable::RHO);
    exponentIdx = (maxResidualOld_[rhoIdx] > maxResidual_[rhoIdx]) ? 1:0;
    
    actualCourantNum_ *= 
      pow(
        maxResidualOld_[rhoIdx]/maxResidual_[rhoIdx], 
        config_->getAdaptiveCourantNumExponent()[exponentIdx]
      );
    
    actualCourantNum_ = max(config_->getCourantNum(), actualCourantNum_);
    actualCourantNum_ = min(config_->getAdaptiveCourantNumMax(), actualCourantNum_);
    
  }
  
}

void SolverImplicit::solve () {
  
  ostringstream residualFileName, residualStr;
  ofstream residualFile;
  
  residualFileName << config_->getOutputFolderPath() << "/residual.csv";
  residualFile.open(residualFileName.str());
  residualFile.close();
  
#ifdef _OPENMP
  cout << "Running solver with " << omp_get_max_threads() << " threads..." << endl << endl;
#else
  cout << "Running solver in serial..." << endl << endl;
#endif
  
  cout << "  N: Time step number." << endl;
  cout << "  W: Conserved variable absolute residual." << endl;
  cout << "  IT: Linear solver number of iterations." << endl;
  cout << "  ERR: Linear solver error." << endl << endl;
  cout << "  | N | ";
  for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
    cout << "W" << (iVariable+1) << " | ";
  }
  cout << "IT | ERR |";
  cout << endl;
  
  for (unsigned long iTimeStep=0; iTimeStep<config_->getNumTimeStep(); iTimeStep++) {
    
    cout << "  | " << iTimeStep+1 << " | ";
    
    residualStr << iTimeStep+1 << ",";
    /*Solver->*/appendToResidual_(residualFileName.str(), residualStr.str());
    residualStr.str("");
    residualStr.clear();
    
    timeStep_(iTimeStep, faceLoopVariables_, residualFileName.str());
    
    if ( 
      iTimeStep%config_->getNumTimeStepBetweenSave() == 0 || 
      iTimeStep == config_->getNumTimeStep()-1 ) {
        
      /*Solver->*/writeToFile_(iTimeStep, config_->getNumCurrentSolution()-1);
      
    }
    
  }
  
  cout << endl << "Solver run completed." << endl;
}