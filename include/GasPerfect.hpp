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
#ifndef GAS_PERFECT_HPP
#define GAS_PERFECT_HPP

#include <cmath>

#include "../include/Gas.hpp"

template <typename T>
class GasPerfect_: public Gas_<T> {
  // Perfect gas class for calculating thermodynamic quantities.
  // Template class to allow use with Eigen autodiff variables.
  
  private:
    
    double gasConstant_, ratioSpecificHeat_;
    double heatCapacityVolume_, heatCapacityPressure_;
    
  public:
    
    // Virtual functions
    //------------------
    void initialize(Config *config);
    T getRatioSpecificHeat() const;
    
    // Primitive variables
    T computeDensity(T const *conservedVariable) const;
    T computeVelocity(T const *conservedVariable, const short &iDim) const;
    T computePressure(T const *conservedVariable) const;
    T computeTemperature(T const *conservedVariable) const;
    void computeAllPrimitiveVariable(
      T const *conservedVariable, T *primitiveVariable) const;
    
    // Other variables
    T computeTemperature(const T pressure, const T density) const;
    T computeVelocitySquared(T const *conservedVariable) const;
    T computeTotalInternalEnergy(T const *conservedVariable) const;
    T computeTotalEnthalpy(T const *conservedVariable) const;
    T computeSoundSpeed(const T temperature) const;
    T computeSoundSpeed(T const *conservedVariable) const;
    T computeContravariantVelocityFromConserved(T const *conservedVariable, double const *nonDimNormal) const;
    T computeMachNum(const T velocity, const T temperature) const;
    
    // Conserved variables
    void computeAllConservedVariable(
      T const *primitiveVariable, T *conservedVariable) const;
    
    // Flux variables
    T computeEnergyFlux(
      T const *primitiveVariable, const T contravariantVelocity) const;
    T computeEnergyFlux(
      T const *primitiveVariable, double const *nonDimNormal) const;
    void computeAllFlux(
      T const *primitiveVariable, double const *nonDimNormal, T *flux) const;
    
    // Convert to dimensional variables
    void convertDimensional(
      T const * primitiveVariable, T * dimensionalPrimitiveVariable);
    
};

#include "../src/GasPerfect.tpp"
typedef GasPerfect_<double> GasPerfect; // Instantiate with double data type for normal use.

#endif // GAS_PERFECT_HPP