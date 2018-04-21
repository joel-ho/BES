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

#include "../include/GasPerfect.hpp"

template <typename T>
void GasPerfect_<T>::initialize (Config *config) {
  this->config_ = config;
  gasConstant_ = this->config_->getGasConstant();
  ratioSpecificHeat_ = this->config_->getRatioSpecificHeat();
  heatCapacityVolume_ = gasConstant_/(ratioSpecificHeat_ - 1.0);
  heatCapacityPressure_ = heatCapacityVolume_ + gasConstant_;
}

template <typename T>
T GasPerfect_<T>::getRatioSpecificHeat () const {
  return ratioSpecificHeat_;
}

template <typename T>
T GasPerfect_<T>::computeDensity (T const *conservedVariable) const {
  return conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO)];
}

template <typename T>
T GasPerfect_<T>::computeVelocity (T const *conservedVariable, const short &iDim) const {
  return 
    (conservedVariable[this->config_->mapConservedVelocityComponent(iDim)]/
    conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO)]);
}

template <typename T>
T GasPerfect_<T>::computePressure (T const *conservedVariable) const {
  
  T tmpValue;
  
  tmpValue = 0;
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    tmpValue += 
      conservedVariable[this->config_->mapConservedVelocityComponent(iDim)]*
      conservedVariable[this->config_->mapConservedVelocityComponent(iDim)];
  }
  
  return
    (ratioSpecificHeat_ - 1.0)*
    (conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO_E)] - 
      tmpValue/
      (2.0*conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO)]));
}

template <typename T>
T GasPerfect_<T>::computeTemperature (T const *conservedVariable) const {
  
  T tmpValue;
  
  tmpValue = 0;
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    tmpValue += 
      conservedVariable[this->config_->mapConservedVelocityComponent(iDim)]*
      conservedVariable[this->config_->mapConservedVelocityComponent(iDim)];
  }
  
  return (
    (1.0/
    (
      heatCapacityVolume_*conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO)]))*
    (conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO_E)] - 
      tmpValue/(2.0*conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO)]))
  );
}

template <typename T>
void GasPerfect_<T>::computeAllPrimitiveVariable (
  T const *conservedVariable, T *primitiveVariable) const {
  
  switch (this->config_->getEquationSet()) {
  
    case EquationSet::EULER_NON_DIMENSIONAL:
      primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::RHO)] = 
        computeDensity(conservedVariable);
      for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
        primitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)] = 
          computeVelocity(conservedVariable, iDim);
      }
      primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::P)] = 
        computePressure(conservedVariable);
      primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::T)] = 
        computeTemperature(conservedVariable);
      break;
      
    default:
      break;
      
  }
  
}

template <typename T>
T GasPerfect_<T>::computeTemperature(const T pressure, const T density) const {
  return pressure/(gasConstant_*density);
}

template <typename T>
T GasPerfect_<T>::computeVelocitySquared (T const *conservedVariable) const {
  
  T tmpValue;
  
  tmpValue = 0;
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    tmpValue += 
      computeVelocity(conservedVariable, iDim)*
      computeVelocity(conservedVariable, iDim);
  }
  return tmpValue;
  
}

template <typename T>
T GasPerfect_<T>::computeTotalInternalEnergy (T const *conservedVariable) const {
  return (conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO_E)]/
    conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO)]
  );
}

template <typename T>
T GasPerfect_<T>::computeTotalEnthalpy (T const *conservedVariable) const {
  return (
    computeTotalInternalEnergy(conservedVariable) + 
    computePressure(conservedVariable)/computeDensity(conservedVariable) 
  );
}

template <typename T>
T GasPerfect_<T>::computeSoundSpeed(const T temperature) const {
  return sqrt(ratioSpecificHeat_*gasConstant_*temperature);
}

template <typename T>
T GasPerfect_<T>::computeSoundSpeed (T const *conservedVariable) const {
  return sqrt(
    ratioSpecificHeat_*(ratioSpecificHeat_ - 1)*
    (computeTotalEnthalpy(conservedVariable) - 
      computeVelocitySquared(conservedVariable)/2.0)
  );
}

template <typename T>
T GasPerfect_<T>::computeContravariantVelocityFromConserved (
  T const *conservedVariable, double const *nonDimNormal) const {
    
  T tmpValue;
  
  tmpValue = 0;
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    tmpValue += 
      computeVelocity(conservedVariable, iDim)*nonDimNormal[iDim];
  }
  return tmpValue;
}

template <typename T>
T GasPerfect_<T>::computeMachNum(const T velocity, const T temperature) const {
  return velocity/computeSoundSpeed(temperature);
}

template <typename T>
void GasPerfect_<T>::computeAllConservedVariable (
  T const *primitiveVariable, T *conservedVariable) const {
    
  T tmpValue;
  tmpValue = 0;
  
  switch (this->config_->getEquationSet()) {
  
    case EquationSet::EULER_NON_DIMENSIONAL:
    
    conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO)] = 
      primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::RHO)];
    
    for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
      conservedVariable[this->config_->mapConservedVelocityComponent(iDim)] =
        conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO)]*
        primitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)];
      tmpValue += 
        primitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)]*
        primitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)];
    }
    
    conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO_E)] = (
      conservedVariable[this->config_->mapConservedVariable(FlowVariable::RHO)]*
      (
        heatCapacityVolume_*primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::T)] 
        + 
        tmpValue/2.0
      )
    );
    
    break;
  }
  
}

template <typename T>
T GasPerfect_<T>::computeEnergyFlux(
  T const *primitiveVariable, const T contravariantVelocity) const {
  
  T velocitySquared = 0;
  
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    velocitySquared += 
      primitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)]*
      primitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)];
  }
  
  return (
    primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::RHO)]*
    (
      heatCapacityPressure_*primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::T)] + 
      velocitySquared/2.0
    )
    *
    contravariantVelocity
  );
  
}

template <typename T>
T GasPerfect_<T>::computeEnergyFlux(T const *primitiveVariable, double const *nonDimNormal) const {
  
  T velocitySquared = 0;
  
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    velocitySquared += 
      primitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)]*
      primitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)];
  }
  
  return (
    primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::RHO)]*
    (
      heatCapacityPressure_*primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::T)] + 
      velocitySquared/2.0
    )
    *
    this->computeContravariantVelocityFromPrimitive(primitiveVariable, nonDimNormal)
  );
  
}

template <typename T>
void GasPerfect_<T>::computeAllFlux(
  T const *primitiveVariable, double const *nonDimNormal, T *flux) const {
    
  T contravariantVelocity;
  
  switch (this->config_->getEquationSet()) {
  
    case EquationSet::EULER_NON_DIMENSIONAL:
    
    contravariantVelocity = this->computeContravariantVelocityFromPrimitive(
      primitiveVariable, nonDimNormal);
    
    flux[this->config_->mapConservedVariable(FlowVariable::RHO)] = 
      this->computeDensityFlux(primitiveVariable, contravariantVelocity);
    
    for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
      flux[this->config_->mapConservedVelocityComponent(iDim)] =
        this->computeMomentumFlux(primitiveVariable, nonDimNormal, contravariantVelocity, iDim);
    }
    flux[this->config_->mapConservedVariable(FlowVariable::RHO_E)] =
      computeEnergyFlux(primitiveVariable, contravariantVelocity);
    
    break;
  }
  
}

template <typename T>
void GasPerfect_<T>::convertDimensional(T const * primitiveVariable, T * dimensionalPrimitiveVariable) {
  
  T soundSpeed;
  
  switch (this->config_->getEquationSet()) {
  
    case EquationSet::EULER_NON_DIMENSIONAL:
      
      dimensionalPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::RHO)] = 
        primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::RHO)];
        
      dimensionalPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::P)] = 
        ratioSpecificHeat_*primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::P)];
        
      dimensionalPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::T)] = 
        ratioSpecificHeat_*
        gasConstant_*
        primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::T)];
        
      soundSpeed = computeSoundSpeed(
        primitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::T)]);
        
      for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
        dimensionalPrimitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)] = 
          primitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)]/soundSpeed;
      }
      
      break;
      
    default:
      break;
  
  }
  
}
