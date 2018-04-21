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

#include "../include/Gas.hpp"

template <typename T>
Gas_<T>::Gas_ () {
  config_ = NULL;
}

template <typename T>
Gas_<T>::~Gas_ () {}

template <typename T>
T Gas_<T>::computeContravariantVelocityFromPrimitive(
  T const *primitiveVariable, double const *nonDimNormal) const {
  
  T tmpValue;
  tmpValue = 0;
  for (short iDim=0; iDim<config_->getNumDim(); iDim++) {
    tmpValue += 
      primitiveVariable[config_->mapPrimitiveVelocityComponent(iDim)]*
      nonDimNormal[iDim];
  }
  return tmpValue;
}

template <typename T>
T Gas_<T>::computeDensityFlux (
  T const *primitiveVariable, const T contravariantVelocity) const {
  return (
    primitiveVariable[config_->mapPrimitiveVariable(FlowVariable::RHO)]*
    contravariantVelocity
  );
}

template <typename T>
T Gas_<T>::computeDensityFlux (T const *primitiveVariable, double const *nonDimNormal) const {
  return (
    primitiveVariable[config_->mapPrimitiveVariable(FlowVariable::RHO)]*
    computeContravariantVelocityFromPrimitive(primitiveVariable, nonDimNormal)
  );
}

template <typename T>
T Gas_<T>::computeMomentumFlux(
  T const *primitiveVariable, double const *nonDimNormal, 
  const T contravariantVelocity, const short &iDim) const {
    
  return (
    primitiveVariable[config_->mapPrimitiveVariable(FlowVariable::RHO)]*
    primitiveVariable[config_->mapPrimitiveVelocityComponent(iDim)]*
    contravariantVelocity
    +
    primitiveVariable[config_->mapPrimitiveVariable(FlowVariable::P)]*
    nonDimNormal[iDim]
  );
  
}

template <typename T>
T Gas_<T>::computeMomentumFlux(
  T const *primitiveVariable, double const *nonDimNormal, const short &iDim) const {
  return (
    primitiveVariable[config_->mapPrimitiveVariable(FlowVariable::RHO)]*
    primitiveVariable[config_->mapPrimitiveVelocityComponent(iDim)]*
    computeContravariantVelocityFromPrimitive(primitiveVariable, nonDimNormal)
    +
    primitiveVariable[config_->mapPrimitiveVariable(FlowVariable::P)]*
    nonDimNormal[iDim]
  );
}
