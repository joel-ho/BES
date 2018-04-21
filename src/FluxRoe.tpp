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

#include "../include/FluxRoe.hpp"

template<typename T>
void FluxRoe_<T>::initialize (Config *config, Gas_<T> *gas) {
  
  this->initializeBase_(config, gas);
  
  nRoeVariable_ = VAR_BUFFER;
  
  // Set up Roe variable map
  roeVariableMap_[ static_cast<short> (FlowVariable::RHO) ] = 
    this->config_->mapPrimitiveVariable(FlowVariable::RHO);
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    roeVariableMap_[ static_cast<short> (this->config_->getPrimitiveVariableVelocityComponent(iDim)) ] = 
      this->config_->mapPrimitiveVelocityComponent(iDim);
  }
  roeVariableMap_[ static_cast<short> (FlowVariable::P) ] = this->config_->mapPrimitiveVariable(FlowVariable::P);
  roeVariableMap_[ static_cast<short> (FlowVariable::T) ] = this->config_->mapPrimitiveVariable(FlowVariable::T);
  
  roeVariableMap_[ static_cast<short> (FlowVariable::H) ] = this->config_->getNumPrimitiveVariable();
  roeVariableMap_[ static_cast<short> (FlowVariable::C) ] = this->config_->getNumPrimitiveVariable() + 1;
  roeVariableMap_[ static_cast<short> (FlowVariable::V_CONTRAVARIANT) ] = this->config_->getNumPrimitiveVariable() + 2;
  roeVariableMap_[ static_cast<short> (FlowVariable::Q_SQ) ] = this->config_->getNumPrimitiveVariable() + 3;
  
}

template<typename T>
FluxRoe_<T>::FluxRoe_() {}

template<typename T>
FluxRoe_<T>::~FluxRoe_ () {}

template<typename T>
short FluxRoe_<T>::mapRoeVariable_(const FlowVariable flowVariable) {
  return roeVariableMap_[ static_cast<short>(flowVariable) ];
}

template<typename T>
short FluxRoe_<T>::mapRoeVelocityComponent_(const short &iDim) {
  return roeVariableMap_[ static_cast<short>(this->config_->getPrimitiveVariableVelocityComponent(iDim)) ];
}

template<typename T>
void FluxRoe_<T>::computePrimitiveVariable_(
  T const *conservedVariable, double const *nonDimNormal, T *primitiveVariable) {
  
  this->gas_->computeAllPrimitiveVariable(conservedVariable, primitiveVariable);
  
  primitiveVariable[mapRoeVariable_(FlowVariable::H)] = 
    this->gas_->computeTotalEnthalpy(conservedVariable);
  primitiveVariable[mapRoeVariable_(FlowVariable::C)] = 
    this->gas_->computeSoundSpeed(conservedVariable);
  primitiveVariable[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)] = 
    this->gas_->computeContravariantVelocityFromConserved(conservedVariable, nonDimNormal);
  primitiveVariable[mapRoeVariable_(FlowVariable::Q_SQ)] = 
    this->gas_->computeVelocitySquared(conservedVariable);
  
}

template<typename T>
T FluxRoe_<T>::computeNumerator_ (const FlowVariable flowVariable) {
  return (
    leftPrimitiveVariable_[mapRoeVariable_(flowVariable)]*
    sqrt(leftPrimitiveVariable_[mapRoeVariable_(FlowVariable::RHO)])
    + 
    rightPrimitiveVariable_[mapRoeVariable_(flowVariable)]*
    sqrt(rightPrimitiveVariable_[mapRoeVariable_(FlowVariable::RHO)])
  );
}

template<typename T>
T FluxRoe_<T>::computeRoeAveragedVelocity_ (
  const short iDim, const T denominator, double const *nonDimNormal) {
  
  roeAveragedVariable_[mapRoeVelocityComponent_(iDim)] = (
    computeNumerator_(this->config_->getPrimitiveVariableVelocityComponent(iDim))
    /
    denominator
  );
  roeAveragedVariable_[mapRoeVariable_(FlowVariable::Q_SQ)] += (
    roeAveragedVariable_[mapRoeVelocityComponent_(iDim)]*
    roeAveragedVariable_[mapRoeVelocityComponent_(iDim)]
  );
  roeAveragedVariable_[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)] += (
    roeAveragedVariable_[mapRoeVelocityComponent_(iDim)]*
    nonDimNormal[iDim]
  );
  
}

template<typename T>
void FluxRoe_<T>::computeRoeAveragedVariable_(double const *nonDimNormal) {
  
  T denominator;

  denominator = (
    sqrt(leftPrimitiveVariable_[mapRoeVariable_(FlowVariable::RHO)]) + 
    sqrt(rightPrimitiveVariable_[mapRoeVariable_(FlowVariable::RHO)])
  );
  
  roeAveragedVariable_[mapRoeVariable_(FlowVariable::Q_SQ)] = 0;
  roeAveragedVariable_[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)] = 0;
  
  roeAveragedVariable_[mapRoeVariable_(FlowVariable::RHO)] = sqrt(
    leftPrimitiveVariable_[mapRoeVariable_(FlowVariable::RHO)]*
    rightPrimitiveVariable_[mapRoeVariable_(FlowVariable::RHO)]
  );
  
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    computeRoeAveragedVelocity_(iDim, denominator, nonDimNormal);
  }
  
  roeAveragedVariable_[mapRoeVariable_(FlowVariable::H)] = (
    computeNumerator_(FlowVariable::H)/denominator
  );
  
  roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)] = sqrt(
    ( this->gas_->getRatioSpecificHeat() - 1 )*
    ( 
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::H)] - 
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::Q_SQ)]/2.0
    )
  );
  
}

template<typename T>
void FluxRoe_<T>::addConvectiveFlux_(
  T const *primitiveVariable, double const *nonDimNormal, T *flux) {
  
  this->gas_->computeAllFlux(primitiveVariable, nonDimNormal, tmpFlux_);
  
  for (short iVariable=0; iVariable<this->config_->getNumConservedVariable(); iVariable++) {
    flux[iVariable] += tmpFlux_[iVariable];
  }
  
}

template<typename T>
T FluxRoe_<T>::hartenCorrection_ (const T eigenvalue) {
  T delta;
  delta = abs(
    roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]*
    this->config_->getHartenCorrectionMachNum()
  );
  if (eigenvalue > delta) {
    return eigenvalue;
  }
  else {
    return (eigenvalue*eigenvalue + delta*delta)/(2.0*delta);
  }
  
}

template<typename T>
void FluxRoe_<T>::subtractDeltaFone_ (double const *nonDimNormal, T *flux) {
  T constant;
  
  constant = (
    hartenCorrection_(abs(
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)] - 
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)] ))
    *
    (
      deltaPrimitiveVariable_[mapRoeVariable_(FlowVariable::P)] - 
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::RHO)]*
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]*
      deltaPrimitiveVariable_[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)]
    )
    /
    (
      2.0*
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]*
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]
    )
  );
  
  flux[this->config_->mapConservedVariable(FlowVariable::RHO)] -= constant;
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    flux[this->config_->mapConservedVelocityComponent(iDim)] -= constant*(
      roeAveragedVariable_[mapRoeVelocityComponent_(iDim)] - 
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]*
      nonDimNormal[iDim]
    );
  }
  flux[this->config_->mapConservedVariable(FlowVariable::RHO_E)] -= constant*(
    roeAveragedVariable_[mapRoeVariable_(FlowVariable::H)] - 
    roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]*
    roeAveragedVariable_[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)]
  );
}

template<typename T>
void FluxRoe_<T>::subtractDeltaFtwo_ (double const *nonDimNormal, T *flux) {
  
  T correctedEigenvalue, constantOne, constantTwo, tmpValue;
  
  tmpValue = 0;
  
  correctedEigenvalue = hartenCorrection_(abs(
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)] 
  ));
  constantOne = (
    correctedEigenvalue
    *
    (deltaPrimitiveVariable_[mapRoeVariable_(FlowVariable::RHO)] - 
      deltaPrimitiveVariable_[mapRoeVariable_(FlowVariable::P)]/
      (
        roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]*
        roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]
      )
    )
  );  
  constantTwo = correctedEigenvalue*
    roeAveragedVariable_[mapRoeVariable_(FlowVariable::RHO)];
  
  flux[this->config_->mapConservedVariable(FlowVariable::RHO)] -= constantOne;
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    
    flux[this->config_->mapConservedVelocityComponent(iDim)] -= (
      constantOne*roeAveragedVariable_[mapRoeVelocityComponent_(iDim)] 
      + 
      constantTwo*
      (
        deltaPrimitiveVariable_[mapRoeVelocityComponent_(iDim)] 
        -
        deltaPrimitiveVariable_[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)]*
        nonDimNormal[iDim]
      )
    );
    
    tmpValue += (
      deltaPrimitiveVariable_[mapRoeVelocityComponent_(iDim)] *
      roeAveragedVariable_[mapRoeVelocityComponent_(iDim)] 
    );
    
  }
  flux[this->config_->mapConservedVariable(FlowVariable::RHO_E)] -= (
    constantOne*roeAveragedVariable_[mapRoeVariable_(FlowVariable::Q_SQ)]/2.0 
    +
    constantTwo*
    (
      tmpValue - 
      deltaPrimitiveVariable_[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)]*
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)]
    )
  );
}

template<typename T>
void FluxRoe_<T>::subtractDeltaFfive_ (double const *nonDimNormal, T *flux) {
  T constant;
  
  constant = (
    hartenCorrection_(abs(
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)] + 
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)] ))
    *
    (
      deltaPrimitiveVariable_[mapRoeVariable_(FlowVariable::P)] + 
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::RHO)]*
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]*
      deltaPrimitiveVariable_[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)]
    )
    /
    (
      2.0*
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]*
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]
    )
  );
  
  flux[this->config_->mapConservedVariable(FlowVariable::RHO)] -= constant;
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    flux[this->config_->mapConservedVelocityComponent(iDim)] -= constant*(
      roeAveragedVariable_[mapRoeVelocityComponent_(iDim)] + 
      roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]*
      nonDimNormal[iDim]
    );
  }
  flux[this->config_->mapConservedVariable(FlowVariable::RHO_E)] -= constant*(
    roeAveragedVariable_[mapRoeVariable_(FlowVariable::H)] + 
    roeAveragedVariable_[mapRoeVariable_(FlowVariable::C)]*
    roeAveragedVariable_[mapRoeVariable_(FlowVariable::V_CONTRAVARIANT)]
  );
}

template<typename T>
void FluxRoe_<T>::computeFlux (
  T const *leftVariable, T const *rightVariable, 
  double const *nonDimNormal, T *flux) {
  // leftVariable and rightVariable are conserved variables
  // flux array is overridden
  
  computePrimitiveVariable_(leftVariable, nonDimNormal, leftPrimitiveVariable_);
  computePrimitiveVariable_(rightVariable, nonDimNormal, rightPrimitiveVariable_);
  for (short iVar=0; iVar<nRoeVariable_; iVar++) {
    deltaPrimitiveVariable_[iVar] = 
      rightPrimitiveVariable_[iVar] - leftPrimitiveVariable_[iVar];
  }
  computeRoeAveragedVariable_(nonDimNormal);
  
  for (short iVar=0; iVar<this->config_->getNumConservedVariable(); iVar++) {
    flux[iVar] = 0;
  }
  addConvectiveFlux_(leftPrimitiveVariable_, nonDimNormal, flux);
  addConvectiveFlux_(rightPrimitiveVariable_, nonDimNormal, flux);
  subtractDeltaFone_(nonDimNormal, flux);
  subtractDeltaFtwo_(nonDimNormal, flux);
  subtractDeltaFfive_(nonDimNormal, flux);
  for (short iVar=0; iVar<this->config_->getNumConservedVariable(); iVar++) {
    flux[iVar] /= 2.0;
  }
  
}
