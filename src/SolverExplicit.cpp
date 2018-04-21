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

#include "../include/SolverExplicit.hpp"

PrivateVariableWrapper::PrivateVariableWrapper() {}

PrivateVariableWrapper::PrivateVariableWrapper(const PrivateVariableWrapper &wrapper) {
  
  config = wrapper.config;
  gas = wrapper.gas;
  nVariable = wrapper.nVariable;
  
  iTimeStep = wrapper.iTimeStep;
  
  fluxCalculator.clear();
  boundaryVariable.clear();
  
  leftVariable = new double [nVariable];
  rightVariable = new double [nVariable];
  flux = new double [nVariable];
  for (short iVariable=0; iVariable<nVariable; iVariable++) {
    leftVariable[iVariable] = wrapper.leftVariable[iVariable];
    rightVariable[iVariable] = wrapper.rightVariable[iVariable];
    flux[iVariable] = wrapper.flux[iVariable];
  }
  
  for (
    map <BoundaryType, Flux*>::const_iterator it=wrapper.fluxCalculator.begin(); 
    it!=wrapper.fluxCalculator.end(); 
    it++) {
    
    switch (it->first) {
    
      case BoundaryType::INTERNAL:
        switch (config->getFluxDiscretizationScheme()) {
          case FluxDiscretizationScheme::ROE:
            fluxCalculator[BoundaryType::INTERNAL] = new FluxRoe;
            fluxCalculator[BoundaryType::INTERNAL]->initialize(config, gas);
            break;
        } 
        boundaryVariable[BoundaryType::INTERNAL] = rightVariable;
        break; 
        
      case BoundaryType::FREESTREAM:
        fluxCalculator[BoundaryType::FREESTREAM] = new FluxFreestream;
        fluxCalculator[BoundaryType::FREESTREAM]->initialize(config, gas);
        boundaryVariable[BoundaryType::FREESTREAM] = config->getFreestreamPrimitiveVariables(iTimeStep);
        break; 
    
      case BoundaryType::SLIP_WALL:
        fluxCalculator[BoundaryType::SLIP_WALL] = new FluxSlipWall;
        fluxCalculator[BoundaryType::SLIP_WALL]->initialize(config, gas);
        boundaryVariable[BoundaryType::SLIP_WALL] = NULL;
        break; 
    
    } 
    
  } // for it=wrapper.fluxCalculator.begin()
  
}

PrivateVariableWrapper::~PrivateVariableWrapper() {
  delete [] leftVariable;
  delete [] rightVariable;
  delete [] flux;
  for (
    map<BoundaryType, Flux*>::iterator it=fluxCalculator.begin(); 
    it!=fluxCalculator.end(); 
    it++) {
    delete it->second;
  }
}


void SolverExplicit::initialize (Config *config, Gas *gas, Mesh *mesh) {
  
  cout << "Initializing explicit solver... " << endl;
  
  BoundaryType activeBoundaryType;
  
  // Save variables
  initializeBase_(config, gas, mesh);
  
  // Save in PrivateVariableWrapper class
  privateVariable_.config = config;
  privateVariable_.gas = gas;
  
  // Initialize storage for internal variables and flux
  privateVariable_.nVariable = config_->getNumConservedVariable();
  privateVariable_.leftVariable = new double [config_->getNumConservedVariable()];
  privateVariable_.rightVariable = new double [config_->getNumConservedVariable()];
  privateVariable_.flux = new double [config_->getNumConservedVariable()];
  for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
    privateVariable_.leftVariable[iVariable] = 0;
    privateVariable_.rightVariable[iVariable] = 0;
    privateVariable_.flux[iVariable] = 0;
  }
  
  // Initialize storage for calculating deltaW
  subResidual_ = new double* [mesh_->getNumElement()];
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    subResidual_[iElement] = new double [config_->getNumConservedVariable()];
    for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
      subResidual_[iElement][iVariable] = 0;
    }
  }
  
  // Initialize storage for maximum residual
  maxResidual_ = new double [config_->getNumConservedVariable()];
  for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
    maxResidual_[iVariable] = 0;
  }
  
  // Initialize required boundary flux calculators
  for (short iBoundary=0; iBoundary<mesh_->getNumActiveBoundary(); iBoundary++) {
    activeBoundaryType = mesh_->getActiveBoundaryType(iBoundary);
    switch (activeBoundaryType) {  
      
      case BoundaryType::FREESTREAM:
        cout << "  Initializing freestream boundary... ";
        privateVariable_.fluxCalculator[activeBoundaryType] = new FluxFreestream;
        privateVariable_.fluxCalculator[activeBoundaryType]->initialize(config_, gas_);
        privateVariable_.boundaryVariable[activeBoundaryType] = 
          config_->getFreestreamPrimitiveVariables(0);
        cout << "Done." << endl;
        break;
        
      case BoundaryType::SLIP_WALL:
        cout << "  Initializing slip wall boundary... ";
        privateVariable_.fluxCalculator[activeBoundaryType] = new FluxSlipWall;
        privateVariable_.fluxCalculator[activeBoundaryType]->initialize(config_, gas_);
        privateVariable_.boundaryVariable[activeBoundaryType] = NULL;
        cout << "Done." << endl;
        break;
    
    }
  }
  
  // Initialize internal flux calculator
  switch (config_->getFluxDiscretizationScheme()) {
  
    case FluxDiscretizationScheme::ROE:
      cout << "  Initializing Roe flux difference splitting... ";
      
      privateVariable_.fluxCalculator[BoundaryType::INTERNAL] = new FluxRoe;
      privateVariable_.fluxCalculator[BoundaryType::INTERNAL]->initialize(config_, gas_);
      break;
      
  }    
  privateVariable_.boundaryVariable[BoundaryType::INTERNAL] = privateVariable_.rightVariable;
  cout << "Done." << endl;
  
  // Initialize gradient calculator
  switch (config_->getGradientScheme()) {
  
    case GradientScheme::GREEN_GAUSS:
      cout << "  Initializing Green-Gauss gradient... ";
      gradient_ = new GradientGreenGauss;
      gradient_->initialize(config_, mesh_, VariableType::CONSERVED);
      cout << "Done." << endl;
      break;
  
  }
  
  // Initialize limiter calculator
  switch (config_->getLimiterScheme()) {
  
    case LimiterScheme::VENKATAKRISHNAN:
      cout << "  Initializing Venkatakrishnan limiter... ";
      limiter_ = new LimiterVenkatakrishnan;
      limiter_->initialize(config_, mesh_, VariableType::CONSERVED);
      cout << "Done." << endl;
      break;
  
  }
  
  // Initialize time step calculator
  switch (config_->getLocalTimeStepOption()) {
  
    case LocalTimeStep::YES:
      cout << "  Initializing local time stepper ";
      timeStepCalculator_ = new TimeStepCalculatorLocal;
      timeStepCalculator_->initialize(config_, gas_, mesh_);
      break;
      
    case LocalTimeStep::NO:
      cout << "  Initializing global time stepper ";
      timeStepCalculator_ = new TimeStepCalculatorGlobal;
      timeStepCalculator_->initialize(config_, gas_, mesh_);
      break;
  
  }
  cout << "with Courant number: " << config_->getCourantNum() << "... ";
  cout << "Done." << endl;
  
  cout << "Solver initialized." << endl << endl;
  
}

SolverExplicit::~SolverExplicit () {
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    delete [] subResidual_[iElement];
  }
  delete [] subResidual_;
  delete [] maxResidual_;
}

void SolverExplicit::clearDeltaConservedVariable_() {
  #pragma omp parallel for default(shared)
  for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
    for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
      subResidual_[iElement][iVariable] = 0;
    }
  }
}

void SolverExplicit::setUniversalInitialCondition(double *primitiveVariable) {
  
  cout << "Initializing flow variables in entire domain... " << endl;
  
  double *conservedVariable;
  
  conservedVariable = new double [config_->getNumConservedVariable()];
  for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
    conservedVariable[iVariable] = 0;
  }
  
  gas_->computeAllConservedVariable(primitiveVariable, conservedVariable);
  
  for (short iSol=0; iSol<config_->getNumCurrentSolution(); iSol++) {
    #pragma omp parallel for default(shared)
    for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      solution_->setConservedVariable(0, iSol, iElement, conservedVariable);
    }
  }
  
  delete [] conservedVariable;
  
  cout << "Flow variables initialized." << endl << endl;
  
}

void SolverExplicit::timeStep_ (
  const unsigned long &iTimeStepStart, PrivateVariableWrapper privateVariable, string residualFileName) {
  
  double newConservedVariable;
  
  SolutionSnapshot *solutionSnapshot;
  
  ostringstream residualStr;
  ofstream residualFile;
  
  timeStepCalculator_->computeAllTimeStep(
    solution_->getCurrentSnapshot(iTimeStepStart, 0),
    config_->getCourantNum());
  
  privateVariable.iTimeStep = iTimeStepStart;
  
  for (short iSub=0; iSub<config_->getNumSubTimeStep(); iSub++) {
    
    clearDeltaConservedVariable_();
    
    solutionSnapshot = solution_->getCurrentSnapshot(iTimeStepStart, min( (short) 1, iSub));
    
    gradient_->computeGradient(solutionSnapshot);
    
    limiter_->computeLimiterTerm(solutionSnapshot, gradient_);
    
    #pragma omp parallel for default(shared) firstprivate(privateVariable)
    for (unsigned long iFace=0; iFace<mesh_->getNumFace(); iFace++) {
      
      privateVariable.face = mesh_->getFace(iFace);
      
      // Extrapolate variables to face, output to leftVariable and rightVariable
      computeFaceConservedVariable_(
        solutionSnapshot, 
        gradient_, 
        limiter_, 
        iFace, 
        privateVariable.leftVariable, 
        privateVariable.rightVariable);
      
      // Compute flux across face, output to flux
      privateVariable.fluxCalculator[privateVariable.face->getBoundaryType()]->computeFlux(
        privateVariable.leftVariable, 
        privateVariable.boundaryVariable[privateVariable.face->getBoundaryType()], 
        privateVariable.face->getNonDimNormal(), 
        privateVariable.flux
      );
      
      // Compute R term where W_new = W_old - (alpha*deltaT/volume) * R
      privateVariable.element = privateVariable.face->getOwnerElement();
      for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
        
        // Check for divergence
        if (privateVariable.flux[iVariable]!=privateVariable.flux[iVariable]) {
          #pragma omp critical 
          {
          cout << endl << "ERROR: Nan values encountered in flux across face at coordinates [";
          for (short iDimError=0; iDimError<mesh_->getNumDim()-1; iDimError++) {
            cout << scientific << privateVariable.face->getCenter()->getSingleCoord(iDimError) << ", ";
          }
          cout << scientific << privateVariable.face->getCenter()->getSingleCoord(mesh_->getNumDim()-1);
          cout << "]. Exiting..." << endl;
          exit(EXIT_FAILURE);
          }
        }
        
        // Update element flux
        #pragma omp atomic
        subResidual_[privateVariable.element->getId()][iVariable] += 
          privateVariable.flux[iVariable]*privateVariable.face->getArea();
        
      }
      
      privateVariable.element = privateVariable.face->getNeighborElement();
      if (privateVariable.element != NULL) {
        for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
          
          // Update element flux
          #pragma omp atomic
          subResidual_[privateVariable.element->getId()][iVariable] -= 
            privateVariable.flux[iVariable]*privateVariable.face->getArea();
            
        }
      }
      
    } // #pragma omp parallel for iFace
    
    if (iSub == config_->getNumSubTimeStep()-1) {
      for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
        maxResidual_[iVariable] = 0;
      }
    }
    
    #pragma omp parallel for default(shared) private(newConservedVariable)
    for (unsigned long iElement=0; iElement<mesh_->getNumElement(); iElement++) {
      for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
        
        if (iSub == config_->getNumSubTimeStep()-1) {
          
          #pragma omp critical
          maxResidual_[iVariable] = 
            (abs(subResidual_[iElement][iVariable])>maxResidual_[iVariable])? 
            abs(subResidual_[iElement][iVariable]):maxResidual_[iVariable];
        }
        
        newConservedVariable = (
          solution_->getConservedVariable(iTimeStepStart, 0, iElement, iVariable) - 
          (
            config_->getTimeStepConstant(iSub)*
            timeStepCalculator_->getTimeStep(iElement)/
            mesh_->getElement(iElement)->getVolume()
          )*
          subResidual_[iElement][iVariable]
        );
        
        solution_->setConservedVariable(
          iTimeStepStart, 1, iElement, iVariable, newConservedVariable);
      
      }
    } // #pragma omp parallel for iElement
    
    if (iSub == config_->getNumSubTimeStep()-1) {
      
      for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
        cout << scientific << maxResidual_[iVariable] << " | ";
        residualStr << scientific << maxResidual_[iVariable] << ",";
      }
      cout << endl;
      residualStr << endl;
      /*Solver->*/appendToResidual_(residualFileName, residualStr.str());
      
    }
    
  } // for iSub
}

void SolverExplicit::solve () {
  
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
  
  cout << "  | Time step | ";
  for (short iVariable=0; iVariable<config_->getNumConservedVariable(); iVariable++) {
    cout << "Conserved variable " << (iVariable+1) << " | ";
  }
  cout << endl;
  
  for (unsigned long iTimeStep=0; iTimeStep<config_->getNumTimeStep(); iTimeStep++) {
    
    cout << "  | " << iTimeStep+1 << " | ";
    
    residualStr << iTimeStep+1 << ",";
    /*Solver->*/appendToResidual_(residualFileName.str(), residualStr.str());
    residualStr.str("");
    residualStr.clear();
    
    timeStep_(iTimeStep, privateVariable_, residualFileName.str());
    
    if ( 
      iTimeStep%config_->getNumTimeStepBetweenSave() == 0 || 
      iTimeStep == config_->getNumTimeStep()-1 ) {
        
      /*Solver->*/writeToFile_(iTimeStep, config_->getNumCurrentSolution()-1);
      
    }
    
  }
  
  cout << endl << "Solver run completed." << endl;
}
