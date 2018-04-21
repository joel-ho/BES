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
#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <iostream>
#include <algorithm>
#include <string>
#include <cmath>

#include "../include/Options.hpp"
#include "../include/FileReader.hpp"

using namespace std;

class Mesh;

class Config {
  
  private:
    
    short nDim_;
    short nConservedVariable_, nPrimitiveVariable_;
    short nCurrentSolution_, nSubTimeStep_;
    short *nVariable_;
    short *conservedMap_, *primitiveMap_;
    
    int linearSolverMaxIter_;
    
    unsigned long nTimeStep_, nTimeStepBetweenSave_, nRampFreestream_;
    
    double gasConstant_, ratioSpecificHeat_;
    double venkatakrishnanLimiterConstant_, hartenCorrectionMachNum_;
    double courantNum_, courantNumMax_, *adaptiveCourantNumExponent_;
    double freestreamRampStartFraction_;
    double *timeStepCoefficient_;
    double *freestreamMachBuffer_, *freestreamPrimitiveVariables_, *rampFreestreamPrimitiveVariables_;
    double linearSolverAbsTol_;
    
    string meshFilePath_, outputFolderPath_;
    
    FlowVariable *conservedVariableVelocityComponents_;
    FlowVariable *primitiveVariableVelocityComponents_;
    
    MeshFormat meshFormat_;
    EquationSet equationSet_;
    SolverType solverType_;
    
    ExplicitScheme explicitScheme_;
    LocalTimeStep localTimeStep_;
    
    GradientScheme gradientScheme_;
    LimiterScheme limiterScheme_;
    FluxDiscretizationScheme fluxDiscretizationScheme_;
    
    AmgPrecondCoarsener amgPrecondCoarsener_;
    AmgPrecondSmoother amgPrecondSmoother_;
    LinearSolver linearSolver_;
    
    OutputFormat outputFormat_;
    OutputSlipWall outputSlipWall_;
    
    void setDefaultValues_();
    void readConfigFile_(string configFilePath);
  
  public:
    
    Config();
    ~Config();
    
    void initialize(string configFilePath);
    void initializePostMesh(Mesh* mesh);
    
    double getGasConstant();
    double getRatioSpecificHeat();
    
    string getMeshFilePath();
    MeshFormat getMeshFormat();
    short getNumDim();
    
    EquationSet getEquationSet();
    short getNumConservedVariable();
    short getNumPrimitiveVariable();
    short getNumVariable(VariableType variableType);
    short mapConservedVariable(FlowVariable flowVariable);
    short mapConservedVelocityComponent(const short &iDim);
    short mapPrimitiveVariable(FlowVariable flowVariable);
    short mapPrimitiveVelocityComponent(const short &iDim);
    short mapVariable(VariableType variableType, FlowVariable flowVariable);
    short mapVelocityComponent(VariableType variableType, const short &iDim);
    
    SolverType getSolverType();
    
    ExplicitScheme getExplicitScheme();
    LocalTimeStep getLocalTimeStepOption();
    unsigned long getNumTimeStep();
    unsigned long getNumTimeStepBetweenSave();
    
    short getNumSubTimeStep();
    double getTimeStepConstant(const short &iStep);
    
    double getCourantNum();
    double* getAdaptiveCourantNumExponent();
    double getAdaptiveCourantNumMax();
    AdaptiveCourantNum getAdaptiveCourantNumOption();
    
    short getNumCurrentSolution();
    
    GradientScheme getGradientScheme();
    LimiterScheme getLimiterScheme();
    FluxDiscretizationScheme getFluxDiscretizationScheme();
    double getVenkatakrishnanLimiterConstant();
    double getHartenCorrectionMachNum();
    
    AmgPrecondCoarsener getAmgPrecondCoarsener();
    AmgPrecondSmoother getAmgPrecondSmoother();
    LinearSolver getLinearSolverType();
    int getLinearSolverNumIter();
    double getLinearSolverTol();
    
    FlowVariable getConservedVariableVelocityComponent(const short &iDim);
    FlowVariable getPrimitiveVariableVelocityComponent(const short &iDim);
    
    double* getFreestreamPrimitiveVariables(const unsigned long &iTimeStep);
    
    OutputFormat getOutputFormat(); 
    string getOutputFolderPath();
    
    OutputSlipWall getOutputSlipWall();
    
};

#endif // CONFIG_HPP