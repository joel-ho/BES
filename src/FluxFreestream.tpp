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

#include "../include/FluxFreestream.hpp"

template <typename T>
void FluxFreestream_<T>::initialize (Config *config, Gas_<T> *gas) {
  
  this->initializeBase_(config, gas);
  
}

template <typename T>
FluxFreestream_<T>::FluxFreestream_() {}

template <typename T>
FluxFreestream_<T>::~FluxFreestream_ () {}

template <typename T>
void FluxFreestream_<T>::computeSupersonicInlet_(
  T const *leftPrimitiveVariable, T const *rightPrimitiveVariable, 
  double const *nonDimNormal, T *flux) {  
  this->gas_->computeAllFlux(rightPrimitiveVariable, nonDimNormal, flux); 
}

template <typename T>
void FluxFreestream_<T>::computeSupersonicOutlet_(
  T const *leftPrimitiveVariable, T const *rightPrimitiveVariable, 
  double const *nonDimNormal, T *flux) {
  this->gas_->computeAllFlux(leftPrimitiveVariable, nonDimNormal, flux); 
}

template <typename T>
void FluxFreestream_<T>::computeSubsonicInlet_(
  T const *leftPrimitiveVariable, T const *rightPrimitiveVariable, 
  double const *nonDimNormal, T *flux) {
  // Linearized Euler equation solved on boundary
  // linearized about interior condition
  
  T leftSoundSpeed, tmpValue, constantOne, constantTwo;
  
  leftSoundSpeed = this->gas_->computeSoundSpeed(
    leftPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::T)]);
    
  // sum( n_i * (u_a - u_d)_i ) = sum( n_i * (u_right - u_left)_i )
  tmpValue = 0;
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    tmpValue += nonDimNormal[iDim]*(
      rightPrimitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)] - 
      leftPrimitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)]
    );
  }
  
  // rho_0 * c_0 = rho_left * c_left
  constantOne = (
    leftPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::RHO)]*
    leftSoundSpeed
  );
  
  // 0.5*( P_a + P_d - rho_0*c_0*(sum(n_i * (u_a - u_d)_i)) ) = 
  // 0.5*( P_right + P_left - constantOne*tmpValue )
  boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::P)] = (
    rightPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::P)] + 
    leftPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::P)] - 
    constantOne*tmpValue
  )/2;
  
  // (P_a - P_b)/(rho_0 * c_0) = (P_right - P_boundary)/constantOne
  constantTwo = (
    (
      rightPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::P)] - 
      boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::P)]
    )
    /
    constantOne
  );
  
  // rho_a + (P_b - P_a)/(c_0*c_0) = rho_right + (P_boundary - P_right)/(c_0*c_0)
  boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::RHO)] = (
    rightPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::RHO)] + 
    (
      boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::P)] - 
      rightPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::P)]
    )
    /
    (
      leftSoundSpeed*leftSoundSpeed
    )
  );
  
  // T = P/(rho*R)
  boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::T)] = (
    this->gas_->computeTemperature(
      boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::P)],
      boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::RHO)]
    )
  );
  
  // (u_b)_i = (u_a)_i - n_i*(P_a - P_b)/(rho_0*c_0) = 
  // (u_boundary)_i = (u_right)_i = n_i*(P_right - P_boundary)/constantOne
  // (u_boundary)_i = (u_right)_i = n_i*constantTwo
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    boundaryPrimitiveVariable_[this->config_->mapPrimitiveVelocityComponent(iDim)] = (
      rightPrimitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)] - 
      constantTwo*nonDimNormal[iDim]
    );
  }
  
  this->gas_->computeAllFlux(boundaryPrimitiveVariable_, nonDimNormal, flux); 
  
}

template <typename T>
void FluxFreestream_<T>::computeSubsonicOutlet_(
  T const *leftPrimitiveVariable, T const *rightPrimitiveVariable, 
  double const *nonDimNormal, T *flux) {
  
  T leftSoundSpeed, constantOne, constantTwo;
  
  leftSoundSpeed = this->gas_->computeSoundSpeed(
    leftPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::T)]);
  
  // rho_0 * c_0 = rho_left * c_left
  constantOne = (
    leftPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::RHO)]*
    leftSoundSpeed
  );
  
  // P_b = P_a => P_boundary = P_right
  boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::P)] = 
    rightPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::P)];
  
  // (P_d - P_b) / (rho_0 * c_0 ) = (P_left - P_boundary)/constantOne
  constantTwo = (
    (
      leftPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::P)] - 
      boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::P)]
    )
    /
    constantOne
  );
  
  // rho_d + (P_b - P_d)/(c_0*c_0) = rho_left + (P_boundary - P_left)/(c_left*c_left)
  boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::RHO)] = (
    leftPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::RHO)] + 
    (
      boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::P)] - 
      leftPrimitiveVariable[this->config_->mapPrimitiveVariable(FlowVariable::P)]
    )
    /
    (
      leftSoundSpeed*leftSoundSpeed
    )
  );
  
  boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::T)] = (
    this->gas_->computeTemperature(
      boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::P)],
      boundaryPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::RHO)]
    )
  );
  
  // (u_d)_i + n_i*( (P_d - P_b)/(rho_0 * c_0 ) ) = 
  // (u_left)_i + n_i*( (P_left - P_boundary)/(rho_0 * c_0 ) ) = 
  // (u_left)_i + n_i*constantTwo
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    boundaryPrimitiveVariable_[this->config_->mapPrimitiveVelocityComponent(iDim)] = (
      leftPrimitiveVariable[this->config_->mapPrimitiveVelocityComponent(iDim)] + 
      constantTwo*nonDimNormal[iDim]
    );
  }
  
  this->gas_->computeAllFlux(boundaryPrimitiveVariable_, nonDimNormal, flux); 
  
}

template <typename T>
void FluxFreestream_<T>::computeFlux(
  T const *leftVariable, T const *rightVariable, 
  double const *nonDimNormal, T *flux) {
  // leftVariable are interior conserved variables 
  //  (extrapolated to boundary in the case of second order spatial discretization)
  // rightVariable are exterior primitive variables
  // nonDimNormal points outwards
  
  T leftMachNum, leftContravariantVelocity;
  
  this->gas_->computeAllPrimitiveVariable(leftVariable, leftPrimitiveVariable_);
  
  leftContravariantVelocity = this->gas_->computeContravariantVelocityFromPrimitive(
    leftPrimitiveVariable_, nonDimNormal);
  leftMachNum = this->gas_->computeMachNum(
    abs(leftContravariantVelocity),
    leftPrimitiveVariable_[this->config_->mapPrimitiveVariable(FlowVariable::T)]
  );
  
  if (leftMachNum>=1 && leftContravariantVelocity<=0) {
    computeSupersonicInlet_(leftPrimitiveVariable_, rightVariable, nonDimNormal, flux);
  }
  else if (leftMachNum>=1 && leftContravariantVelocity>0) {
    computeSupersonicOutlet_(leftPrimitiveVariable_, rightVariable, nonDimNormal, flux);
  }
  else if (leftMachNum<1 && leftContravariantVelocity<=0) {
    computeSubsonicInlet_(leftPrimitiveVariable_, rightVariable, nonDimNormal, flux);
  }
  else if (leftMachNum<1 && leftContravariantVelocity>0) {
    computeSubsonicOutlet_(leftPrimitiveVariable_, rightVariable, nonDimNormal, flux);
  }
  
  
}
