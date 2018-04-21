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
#ifndef GAS_HPP
#define GAS_HPP

#include "../include/Config.hpp"

template <typename T>
class Gas_ {
  // Gas base class for calculating thermodynamic quantities 
  // (i.e. primitive variables from conserved, flux from 
  // primitive variables etc.). Template to allow use with 
  // Eigen autodiff variables. 
  
  protected:
    
    Config *config_;
    
  public:
  
    Gas_();
    virtual ~Gas_();
    
    // Virtual functions
    //------------------
    virtual void initialize(Config *config) = 0;
    virtual T getRatioSpecificHeat() const = 0;
    
    // Compute primitive variables from conserved variables
    virtual T computeDensity(T const *conservedVariable) const = 0;
    virtual T computeVelocity(
      T const *conservedVariable, const short &iDim) const = 0;
    virtual T computePressure(T const *conservedVariable) const = 0;
    virtual T computeTemperature(T const *conservedVariable) const = 0;
    virtual void computeAllPrimitiveVariable(
      T const *conservedVariable, T *primitiveVariable) const = 0;
    
    // Other variables
    virtual T computeTemperature(const T pressure, const T density) const = 0;
    virtual T computeVelocitySquared(T const *conservedVariable) const = 0;
    virtual T computeTotalInternalEnergy(T const *conservedVariable) const = 0;
    virtual T computeTotalEnthalpy(T const *conservedVariable) const = 0;
    virtual T computeSoundSpeed(const T temperature) const = 0;
    virtual T computeSoundSpeed(T const *conservedVariable) const = 0;
    virtual T computeContravariantVelocityFromConserved(
      T const *conservedVariable, double const *nonDimNormal) const = 0;
    virtual T computeMachNum(
      const T velocity, const T temperature) const = 0;
    
    virtual void computeAllConservedVariable(
      T const *primitiveVariable, T *conservedVariable) const = 0;
    
    virtual T computeEnergyFlux(
      T const *primitiveVariable, const T contravariantVelocity) const = 0;
    virtual T computeEnergyFlux(
      T const *primitiveVariable, double const *nonDimNormal) const = 0;
    virtual void computeAllFlux(
      T const *primitiveVariable, double const *nonDimNormal, T *flux) const = 0;
      
    virtual void convertDimensional(
      T const * primitiveVariable, T * dimensionalPrimitiveVariable) = 0;
    
    // Class methods
    //--------------
    // Compute flux variables from primitive variables
    T computeContravariantVelocityFromPrimitive(
      T const *primitiveVariable, double const *nonDimNormal) const;
    T computeDensityFlux(
      T const *primitiveVariable, const T contravariantVelocity) const;
    T computeDensityFlux(
      T const *primitiveVariable, double const *nonDimNormal) const;
    T computeMomentumFlux(
      T const *primitiveVariable, double const *nonDimNormal, 
      const T contravariantVelocity, const short &iDim) const;
    T computeMomentumFlux(
      T const *primitiveVariable, double const *nonDimNormal, const short &iDim) const;

};

#include "../src/Gas.tpp"
typedef Gas_<double> Gas; // Instantiate with double data type for normal use.

#endif // GAS_HPP