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
#ifndef FLUX_FREESTREAM_HPP
#define FLUX_FREESTREAM_HPP

#include <cmath>
#include <map>

#include "../include/Flux.hpp"
#include "../include/GasPerfect.hpp"

template <typename T>
class FluxFreestream_: public Flux_<T> {
  
  private:
    
    T leftPrimitiveVariable_[VAR_BUFFER], boundaryPrimitiveVariable_[VAR_BUFFER];
    
    // Private functions
    void computeSupersonicInlet_(
      T const *leftPrimitiveVariable, T const *rightPrimitiveVariable, 
      double const *nonDimNormal, T *flux);
    void computeSupersonicOutlet_(
      T const *leftPrimitiveVariable, T const *rightPrimitiveVariable, 
      double const *nonDimNormal, T *flux);
    void computeSubsonicInlet_(
      T const *leftPrimitiveVariable, T const *rightPrimitiveVariable, 
      double const *nonDimNormal, T *flux);
    void computeSubsonicOutlet_(
      T const *leftPrimitiveVariable, T const *rightPrimitiveVariable, 
      double const *nonDimNormal, T *flux);
    
  public:
  
    FluxFreestream_();
    ~FluxFreestream_();
  
    // Virtual functions
    void initialize(Config *config, Gas_<T> *gas);
    void computeFlux(
      T const *leftVariable, T const *rightVariable, 
      double const *nonDimNormal, T *flux);
    
};

#include "../src/FluxFreestream.tpp"
typedef FluxFreestream_<double> FluxFreestream; // Instantiate with T data type for normal use.

#endif // FLUX_FREESTREAM_HPP