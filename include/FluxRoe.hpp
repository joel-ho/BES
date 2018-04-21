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
#ifndef FLUX_ROE_HPP
#define FLUX_ROE_HPP

#include <cmath>

#include "../include/Flux.hpp"
#include "../include/GasPerfect.hpp"

template<typename T>
class FluxRoe_: public Flux_<T> {
  
  private:
    
    // Private variables
    short nRoeVariable_;
    short roeVariableMap_[static_cast<short>(FlowVariable::COUNT)];
    
    T leftPrimitiveVariable_[VAR_BUFFER], rightPrimitiveVariable_[VAR_BUFFER];
    T deltaPrimitiveVariable_[VAR_BUFFER], roeAveragedVariable_[VAR_BUFFER];
    T tmpFlux_[VAR_BUFFER];
    
    // Private functions
    short mapRoeVariable_(const FlowVariable flowVariable);
    short mapRoeVelocityComponent_(const short &iDim);
    void computePrimitiveVariable_(
      T const *conservedVariable, double const *nonDimNormal,
      T *primitiveVariable);
    T computeNumerator_(const FlowVariable flowVariable);
    T computeRoeAveragedVelocity_(
      const short iDim, const T denominator, double const *nonDimNormal);
    void computeRoeAveragedVariable_(double const *nonDimNormal);
    void addConvectiveFlux_(
      T const *primitiveVariable, double const *nonDimNormal, T *flux);
    T hartenCorrection_(const T eigenvalue);
    void subtractDeltaFone_(double const *nonDimNormal, T *flux);
    void subtractDeltaFtwo_(double const *nonDimNormal, T *flux);
    void subtractDeltaFfive_(double const *nonDimNormal, T *flux);
    
  public:
    
    FluxRoe_();
    ~FluxRoe_();
    
    // Virtual functions
    void initialize(Config *config, Gas_<T> *gas);
    void computeFlux(
      T const *leftVariable, T const *rightVariable, 
      double const *nonDimNormal, T *flux);
      
};

#include "../src/FluxRoe.tpp"
typedef FluxRoe_<double> FluxRoe; // Instantiate with T data type for normal use.

#endif // FLUX_ROE_HPP