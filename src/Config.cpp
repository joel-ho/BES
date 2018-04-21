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

#include "../include/Config.hpp"
#include "../include/Mesh.hpp"

Config::Config() {
  
  nDim_ = 0;
  nConservedVariable_ = 0;
  nPrimitiveVariable_ = 0;
  
  nCurrentSolution_ = 0;
  nSubTimeStep_ = 0;
  
  nVariable_ = NULL;
  
  nRampFreestream_ = 0;
  freestreamRampStartFraction_ = 0;
  
  nTimeStep_ = 0;
  nTimeStepBetweenSave_ = 0;
  
  gasConstant_ = 0;
  ratioSpecificHeat_ = 0;
  
  venkatakrishnanLimiterConstant_ = 0;
  hartenCorrectionMachNum_ = 0;
  
  courantNum_ = 0;
  courantNumMax_ = 0;
  adaptiveCourantNumExponent_ = NULL;
  
  timeStepCoefficient_ = NULL;
  freestreamMachBuffer_ = NULL;
  freestreamPrimitiveVariables_ = NULL;
  rampFreestreamPrimitiveVariables_ = NULL;
  
  primitiveMap_ = new short [ static_cast<short> (FlowVariable::COUNT) ];
  conservedMap_ = new short [ static_cast<short> (FlowVariable::COUNT) ];
  for (short iVariable=0; iVariable< static_cast<short>(FlowVariable::COUNT); iVariable++) {
    primitiveMap_[iVariable] = 0;
    conservedMap_[iVariable] = 0;
  }
  
  conservedVariableVelocityComponents_ = new FlowVariable [DIM_BUFFER];
  conservedVariableVelocityComponents_[0] = FlowVariable::RHO_U;
  conservedVariableVelocityComponents_[1] = FlowVariable::RHO_V;
  conservedVariableVelocityComponents_[2] = FlowVariable::RHO_W;
  
  primitiveVariableVelocityComponents_ = new FlowVariable [DIM_BUFFER];
  primitiveVariableVelocityComponents_[0] = FlowVariable::U;
  primitiveVariableVelocityComponents_[1] = FlowVariable::V;
  primitiveVariableVelocityComponents_[2] = FlowVariable::W;
  
}

void Config::setDefaultValues_ () {
  
  outputFormat_ = OutputFormat::VTK;
  
  equationSet_ = EquationSet::EULER_NON_DIMENSIONAL;
  solverType_ = SolverType::EXPLICIT;
  explicitScheme_ = ExplicitScheme::MULTI_STEP;
  localTimeStep_ = LocalTimeStep::NO;
  
  nRampFreestream_ = 0;
  freestreamRampStartFraction_ = 0.5;
  
  courantNum_ = 2.0;
  courantNumMax_ = 0;
  adaptiveCourantNumExponent_ = new double[2];
  adaptiveCourantNumExponent_[0] = 1;
  adaptiveCourantNumExponent_[1] = 1;
  
  nSubTimeStep_ = 4;
  timeStepCoefficient_ = new double [nSubTimeStep_];
  timeStepCoefficient_[0] = 0.0833;
  timeStepCoefficient_[1] = 0.2069;
  timeStepCoefficient_[2] = 0.4265;
  timeStepCoefficient_[3] = 1.0;
  
  gradientScheme_ = GradientScheme::GREEN_GAUSS;
  limiterScheme_ = LimiterScheme::VENKATAKRISHNAN;
  fluxDiscretizationScheme_ = FluxDiscretizationScheme::ROE;
  
  linearSolverMaxIter_ = 1000;
  linearSolverAbsTol_ = 1e-7;
  linearSolver_ = LinearSolver::FGMRES;
  
  venkatakrishnanLimiterConstant_ = 5.0; // Lower number for more limiting - more difficult convergence
  hartenCorrectionMachNum_ = 0.05;
  
  gasConstant_ = 287;
  ratioSpecificHeat_ = 1.4;
  
  outputSlipWall_ = OutputSlipWall::NO;
  
}

void Config::readConfigFile_(string configFilePath) {
  
  double argumentNum[LINE_BUFFER];  
  string fileLine;
  string splitComment[LINE_BUFFER], splitLine[LINE_BUFFER], argument[LINE_BUFFER];
  string commentDelimiter="//";
  string optionDelimiter="=";
  string argumentDelimiter=",";
  
  ifstream configFile;
  
  for (int iEntry=0; iEntry<LINE_BUFFER; iEntry++) {
    splitLine[iEntry] = "";
    argument[iEntry] = "";
    argumentNum[iEntry] = 0;
  }
  
  configFile.open(configFilePath);
  while ( getline(configFile, fileLine) ) {
    
    ReadFile::splitString(fileLine, commentDelimiter, splitComment);
    ReadFile::splitString(splitComment[0], "=", splitLine);
    
    if (splitLine[0].compare("MESH_FORMAT")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("GMSH")==0) {
        meshFormat_ = MeshFormat::GMSH;
      }
      else {
        cout << "ERROR: Mesh format not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("MESH_FILE")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      meshFilePath_ = argument[0];
    }
    
    else if (splitLine[0].compare("OUTPUT_FORMAT")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("VTK")==0) {
        outputFormat_ = OutputFormat::VTK;
      }
      else if (argument[0].compare("CSV")==0) {
        outputFormat_ = OutputFormat::CSV;
      }
      else {
        cout << "ERROR: Output format not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("OUTPUT_SLIP_WALL")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("YES")==0) {
        outputSlipWall_ = OutputSlipWall::YES;
      }
      else if (argument[0].compare("NO")==0) {
        outputSlipWall_ = OutputSlipWall::NO;
      }
      else {
        cout << "ERROR: Output slip wall format not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("OUTPUT_FOLDER")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      outputFolderPath_ = argument[0];
    }
    
    else if (splitLine[0].compare("EQUATION")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("EULER_NON_DIMENSIONAL")==0) {
        equationSet_ = EquationSet::EULER_NON_DIMENSIONAL;
      }
      else {
        cout << "ERROR: Equation not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("SOLVER_TYPE")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("EXPLICIT")==0) {
        solverType_ = SolverType::EXPLICIT;
      }
      else if (argument[0].compare("IMPLICIT")==0) {
        solverType_ = SolverType::IMPLICIT;
      }
      else {
        cout << "ERROR: Solver type not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("EXPLICIT_SCHEME")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("MULTI_STEP")==0) {
        explicitScheme_ = ExplicitScheme::MULTI_STEP;
      }
      else {
        cout << "ERROR: Explicit scheme not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("MULTI_TIME_STEP_OPTIONS")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      nSubTimeStep_ = argumentNum[0];
      delete [] timeStepCoefficient_;
      timeStepCoefficient_ = new double [nSubTimeStep_];
      for (short iCoeff=0; iCoeff<nSubTimeStep_; iCoeff++) {
        timeStepCoefficient_[iCoeff] = argumentNum[iCoeff+1];
      }
    }
    
    else if (splitLine[0].compare("COURANT_NUM")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      courantNum_ = argumentNum[0];
    }
    
    else if (splitLine[0].compare("ADAPTIVE_COURANT_NUM_MAX")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      courantNumMax_ = argumentNum[0];
    }
    
    else if (splitLine[0].compare("ADAPTIVE_COURANT_NUM_EXPONENT")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      adaptiveCourantNumExponent_[0] = argumentNum[0];
      adaptiveCourantNumExponent_[1] = argumentNum[1];
    }
    
    else if (splitLine[0].compare("LOCAL_TIME_STEP")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("YES")==0) {
        localTimeStep_ = LocalTimeStep::YES;
      }
      else if (argument[0].compare("NO")==0) {
        localTimeStep_ = LocalTimeStep::NO;
      }
      else {
        cout << "ERROR: Local time step option not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("NUMBER_OF_TIME_STEP")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      nTimeStep_ = argumentNum[0];
    }
    
    else if (splitLine[0].compare("SAVE_INTERVAL")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      nTimeStepBetweenSave_ = argumentNum[0];
    }
    
    else if (splitLine[0].compare("FREESTREAM_RAMP")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      nRampFreestream_ = argumentNum[0];
    }
    
    else if (splitLine[0].compare("FREESTREAM_RAMP_START_FRAC")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      freestreamRampStartFraction_ = argumentNum[0];
    }
    
    else if (splitLine[0].compare("GRADIENT_SCHEME")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("GREEN_GAUSS")==0) {
        gradientScheme_ = GradientScheme::GREEN_GAUSS;
      }
      else {
        cout << "ERROR: Gradient scheme not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("LIMITER_SCHEME")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("VENKATAKRISHNAN")==0) {
        limiterScheme_ = LimiterScheme::VENKATAKRISHNAN;
      }
      else {
        cout << "ERROR: Limiter not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("VENKATAKRISHNAN_FACTOR")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      venkatakrishnanLimiterConstant_ = argumentNum[0];
    }
    
    else if (splitLine[0].compare("HARTEN_CORRECTION_LIMIT")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      hartenCorrectionMachNum_ = argumentNum[0];
    }
    
    else if (splitLine[0].compare("FLUX_SCHEME")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("ROE")==0) {
        fluxDiscretizationScheme_ = FluxDiscretizationScheme::ROE;
      }
      else {
        cout << "ERROR: Flux discretization scheme not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("GAS_CONSTANT")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      gasConstant_ = argumentNum[0];
    }
    
    else if (splitLine[0].compare("SPECIFIC_HEAT_RATIO")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      ratioSpecificHeat_ = argumentNum[0];
    }
    
    else if (splitLine[0].compare("FREESTREAM_MACH")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      freestreamMachBuffer_ = new double [DIM_BUFFER];
      for (short iDim=0; iDim<DIM_BUFFER; iDim++) {
        freestreamMachBuffer_[iDim] = argumentNum[iDim];
      }
    }
    
    else if (splitLine[0].compare("LINEAR_SOLVER")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("FGMRES")==0) {
        linearSolver_ = LinearSolver::FGMRES;
      }
      else if (argument[0].compare("GMRES")==0) {
        linearSolver_ = LinearSolver::GMRES;
      }
      else if (argument[0].compare("LGMRES")==0) {
        linearSolver_ = LinearSolver::LGMRES;
      }
      else if (argument[0].compare("BICGSTAB")==0) {
        linearSolver_ = LinearSolver::BICGSTAB;
      }
      else {
        cout << "ERROR: Linear solver option not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("AMG_PRECOND_COARSENER")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("AGGREGATION")==0) {
        amgPrecondCoarsener_ = AmgPrecondCoarsener::AGGREGATION;
      }
      else if (argument[0].compare("SMOOTHED_AGGREGATION")==0) {
        amgPrecondCoarsener_ = AmgPrecondCoarsener::SMOOTHED_AGGREGATION;
      }
      else if (argument[0].compare("SMOOTHED_AGGREGATION_EMIN")==0) {
        amgPrecondCoarsener_ = AmgPrecondCoarsener::SMOOTHED_AGGREGATION_EMIN;
      }
      else if (argument[0].compare("RUGE_STUBEN")==0) {
        amgPrecondCoarsener_ = AmgPrecondCoarsener::RUGE_STUBEN;
      }
      else {
        cout << "ERROR: AMG preconditioner coarsening option not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("AMG_PRECOND_SMOOTHER")==0) {
      ReadFile::splitString(splitLine[1], argumentDelimiter, argument);
      if (argument[0].compare("DAMPED_JACOBI")==0) {
        amgPrecondSmoother_ = AmgPrecondSmoother::DAMPED_JACOBI;
      }
      else if (argument[0].compare("ILU0")==0) {
        amgPrecondSmoother_ = AmgPrecondSmoother::ILU0;
      }
      else if (argument[0].compare("SPAI0")==0) {
        amgPrecondSmoother_ = AmgPrecondSmoother::SPAI0;
      }
      else if (argument[0].compare("GAUSS_SEIDEL")==0) {
        amgPrecondSmoother_ = AmgPrecondSmoother::GAUSS_SEIDEL;
      }
      else {
        cout << "ERROR: AMG preconditioner smoothening option not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    else if (splitLine[0].compare("LINEAR_SOLVER_MAX_ITER")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      linearSolverMaxIter_ = argumentNum[0];
    }
    
    else if (splitLine[0].compare("LINEAR_SOLVER_ABS_TOL")==0) {
      ReadFile::splitStringToNum(splitLine[1], argumentDelimiter, argumentNum);
      linearSolverAbsTol_ = argumentNum[0];
    }
    
    else if (splitLine[0].compare("")!=0) {
      cout << "ERROR: Option \""<< splitLine[0] <<"\" not supported. Exiting..." << endl;
      exit(EXIT_FAILURE);
    }
    
  } // while ( getline(configFile, fileLine) )
  
}

void Config::initialize(string configFilePath) {
  
  setDefaultValues_();
  readConfigFile_(configFilePath);
  
  // Get number of current snapshots to store
  if (solverType_ == SolverType::EXPLICIT) {
    
    switch ( explicitScheme_ ) {
      case ExplicitScheme::MULTI_STEP:
        nCurrentSolution_ = 2;
        break;
      default:
        cout << "ERROR: Explicit scheme entered not supported. Exiting..." << endl;
        exit(EXIT_FAILURE);
        break;
    }
    
  }

  else if (solverType_ == SolverType::IMPLICIT) {
    nCurrentSolution_ = 1;
  }
  
  else {
    cout << "ERROR: Solver type entered not supported. Exiting..." << endl;
    exit(EXIT_FAILURE);
  }
  
}

void Config::initializePostMesh (Mesh *mesh) {
  
  nDim_ = mesh->getNumDim();
  
  // Get number of variables and set up variable name mapping
  switch (equationSet_) {
    
    case EquationSet::EULER_NON_DIMENSIONAL:
    
      // Get number of conserved and primitive variables
      nConservedVariable_ = mesh->getNumDim() + 2;
      nPrimitiveVariable_ = mesh->getNumDim() + 3;

      // Set up mapping of flow variable to array index
      conservedMap_[ static_cast<short> (FlowVariable::RHO) ] = 0;
      primitiveMap_[ static_cast<short> (FlowVariable::RHO) ] = 0;
      for (short iDim=0; iDim<mesh->getNumDim(); iDim++) {
        conservedMap_[ static_cast<short> (conservedVariableVelocityComponents_[iDim]) ] = iDim + 1;
        primitiveMap_[ static_cast<short> (primitiveVariableVelocityComponents_[iDim]) ] = iDim + 1;
      }
      conservedMap_[ static_cast<short> (FlowVariable::RHO_E) ] = nConservedVariable_ - 1;
      primitiveMap_[ static_cast<short> (FlowVariable::P) ] = nPrimitiveVariable_ - 2;
      primitiveMap_[ static_cast<short> (FlowVariable::T) ] = nPrimitiveVariable_ - 1;
      
      // Set up freestream primitive variables 
      if (freestreamMachBuffer_!=NULL) {
        freestreamPrimitiveVariables_ = new double [nPrimitiveVariable_];
        freestreamPrimitiveVariables_[0] = 1.0;
        for (short iDim=0; iDim<mesh->getNumDim(); iDim++) {
          freestreamPrimitiveVariables_[iDim+1] = freestreamMachBuffer_[iDim];
        }
        freestreamPrimitiveVariables_[nPrimitiveVariable_ - 2] = 1.0/ratioSpecificHeat_;
        freestreamPrimitiveVariables_[nPrimitiveVariable_ - 1] = 1.0/(ratioSpecificHeat_*gasConstant_);
      }
      
      break;
      
    default:
      cout << "ERROR: Equation set entered not supported. Exiting..." << endl;
      exit(EXIT_FAILURE);
      break;
  
  }
  
  // Create freestream variable ramp
  if (freestreamPrimitiveVariables_ != NULL && nRampFreestream_>0) {
    rampFreestreamPrimitiveVariables_ = new double [nPrimitiveVariable_];
    for (short iVariable=0; iVariable<nPrimitiveVariable_; iVariable++) {
      rampFreestreamPrimitiveVariables_[iVariable] = freestreamPrimitiveVariables_[iVariable];
    }
    for (short iDim=0; iDim<nDim_; iDim++) {
      rampFreestreamPrimitiveVariables_[ mapPrimitiveVelocityComponent(iDim) ] = 
        freestreamRampStartFraction_*freestreamPrimitiveVariables_[ mapPrimitiveVelocityComponent(iDim) ];
    }
  }
  
  nVariable_ = new short [ VariableType::COUNT];
  nVariable_[VariableType::CONSERVED] = nConservedVariable_;
  nVariable_[VariableType::PRIMITIVE] = nPrimitiveVariable_;
  
}

Config::~Config() {
  delete [] nVariable_;
  delete [] conservedMap_;
  delete [] primitiveMap_;
  delete [] conservedVariableVelocityComponents_;
  delete [] primitiveVariableVelocityComponents_;
  if (rampFreestreamPrimitiveVariables_ != NULL) { delete [] rampFreestreamPrimitiveVariables_; }
  if (freestreamMachBuffer_!=NULL) { delete [] freestreamMachBuffer_; }
  if (adaptiveCourantNumExponent_ != NULL) { delete [] adaptiveCourantNumExponent_; }
}

double Config::getGasConstant () {
  return gasConstant_;
}

double Config::getRatioSpecificHeat () {
  return ratioSpecificHeat_;
}

string Config::getMeshFilePath () {
  return meshFilePath_;
}

MeshFormat Config::getMeshFormat () {
  return meshFormat_;
}

short Config::getNumDim () {
  return nDim_;
}

EquationSet Config::getEquationSet() {
  return equationSet_;
}

short Config::getNumConservedVariable () {
  return nConservedVariable_;
}

short Config::getNumPrimitiveVariable () {
  return nPrimitiveVariable_;
}

short Config::getNumVariable(VariableType variableType) {
  return nVariable_[variableType];
}

short Config::mapConservedVariable(FlowVariable flowVariable) {
  return conservedMap_[ static_cast<short> (flowVariable) ];
}

short Config::mapPrimitiveVariable(FlowVariable flowVariable) {
  return primitiveMap_[ static_cast<short> (flowVariable) ];
}

short Config::mapConservedVelocityComponent (const short &iDim) {
  return conservedMap_[ static_cast<short> (conservedVariableVelocityComponents_[iDim]) ];
}

short Config::mapPrimitiveVelocityComponent (const short &iDim) {
  return primitiveMap_[ static_cast<short> (primitiveVariableVelocityComponents_[iDim]) ];
}

short Config::mapVariable(VariableType variableType, FlowVariable flowVariable) {
  return (variableType == VariableType::CONSERVED)? 
    conservedMap_[ static_cast<short> (flowVariable) ]:primitiveMap_[ static_cast<short> (flowVariable) ];
}

short Config::mapVelocityComponent(VariableType variableType, const short &iDim) {
  return (variableType == VariableType::CONSERVED)? 
    conservedMap_[ static_cast<short> (conservedVariableVelocityComponents_[iDim]) ]:
    primitiveMap_[ static_cast<short> (primitiveVariableVelocityComponents_[iDim]) ];
}

SolverType Config::getSolverType() {
  return solverType_;
}

ExplicitScheme Config::getExplicitScheme() {
  return explicitScheme_;
}

LocalTimeStep Config::getLocalTimeStepOption() {
  return localTimeStep_;
}

unsigned long Config::getNumTimeStep() {
  return nTimeStep_;
}

unsigned long Config::getNumTimeStepBetweenSave() {
  return nTimeStepBetweenSave_;
}

short Config::getNumSubTimeStep() {
  return nSubTimeStep_;
}

double Config::getTimeStepConstant(const short &iStep) {
  return timeStepCoefficient_[iStep];
}

double Config::getCourantNum() {
  return courantNum_;
}

double* Config::getAdaptiveCourantNumExponent () {
  return adaptiveCourantNumExponent_;
}

double Config::getAdaptiveCourantNumMax() {
  return courantNumMax_;
}

AdaptiveCourantNum Config::getAdaptiveCourantNumOption() {
  return (courantNumMax_ > courantNum_)? 
    AdaptiveCourantNum::YES : AdaptiveCourantNum::NO;
}

short Config::getNumCurrentSolution () {
  return nCurrentSolution_;
}

GradientScheme Config::getGradientScheme () {
  return gradientScheme_;
}

LimiterScheme Config::getLimiterScheme () {
  return limiterScheme_;
}

FluxDiscretizationScheme Config::getFluxDiscretizationScheme() {
  return fluxDiscretizationScheme_;
}

double Config::getVenkatakrishnanLimiterConstant() {
  return venkatakrishnanLimiterConstant_;
}

double Config::getHartenCorrectionMachNum () {
  return hartenCorrectionMachNum_;
}

AmgPrecondCoarsener Config::getAmgPrecondCoarsener() {
  return amgPrecondCoarsener_;
}

AmgPrecondSmoother Config::getAmgPrecondSmoother() {
  return amgPrecondSmoother_;
}

LinearSolver Config::getLinearSolverType() {
  return linearSolver_;
}

int Config::getLinearSolverNumIter() {
  return linearSolverMaxIter_;
}

double Config::getLinearSolverTol() {
  return linearSolverAbsTol_;
}

FlowVariable Config::getConservedVariableVelocityComponent(const short &iDim) {
  return conservedVariableVelocityComponents_[iDim];
}

FlowVariable Config::getPrimitiveVariableVelocityComponent(const short &iDim) {
  return primitiveVariableVelocityComponents_[iDim];
}

double* Config::getFreestreamPrimitiveVariables(const unsigned long &iTimeStep) {
  
  if (iTimeStep >= nRampFreestream_) {
    return freestreamPrimitiveVariables_;
  }
  else {
    for (short iDim=0; iDim<nDim_; iDim++) {
      rampFreestreamPrimitiveVariables_[ mapPrimitiveVelocityComponent(iDim) ] = 
        freestreamPrimitiveVariables_[ mapPrimitiveVelocityComponent(iDim) ]*
        (
          freestreamRampStartFraction_ + 
          static_cast<double>(iTimeStep)/static_cast<double>(nRampFreestream_)*
            (1-freestreamRampStartFraction_)
        );
    }
    return rampFreestreamPrimitiveVariables_;
  }
  
}

OutputFormat Config::getOutputFormat () {
  return outputFormat_;
}

string Config::getOutputFolderPath() {
  return outputFolderPath_;
}

OutputSlipWall Config::getOutputSlipWall() {
  return outputSlipWall_;
}