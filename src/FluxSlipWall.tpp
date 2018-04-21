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

#include "../include/FluxSlipWall.hpp"

template <typename T>
void FluxSlipWall_<T>::initialize (Config *config, Gas_<T> *gas) {
  this->initializeBase_(config, gas);
}

template <typename T>
FluxSlipWall_<T>::FluxSlipWall_() {}

template <typename T>
FluxSlipWall_<T>::~FluxSlipWall_() {}

template <typename T>
void FluxSlipWall_<T>::computeFlux(
  T const *leftVariable, T const *rightVariable, 
  double const *nonDimNormal, T *flux) {
  // leftVariable are internal conserved variables extrapolated to wall boundary
  // nonDimNormal points inside out of internal element
  // rightVariable can be a NULL pointer
  
  T pressure;
  
  for (short iVariable=0; iVariable<this->config_->getNumConservedVariable(); iVariable++) {
    flux[iVariable] = 0;
  }
  
  pressure = this->gas_->computePressure(leftVariable);
  for (short iDim=0; iDim<this->config_->getNumDim(); iDim++) {
    flux[this->config_->mapConservedVelocityComponent(iDim)] = ( 
      pressure*nonDimNormal[iDim]
    );
  }
  
}
