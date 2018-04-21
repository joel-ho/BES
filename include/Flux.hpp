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
#ifndef FLUX_HPP
#define FLUX_HPP

#include "../include/Config.hpp"
#include "../include/Gas.hpp"

template <typename T>
class Flux_ {
  
  protected:
    
    Config *config_;
    Gas_<T> *gas_;
    
    void initializeBase_(Config *config, Gas_<T> *gas_);
    
  public:
    
    // Constructors and destructors
    Flux_();
    virtual ~Flux_();
    
    // Virtual functions
    virtual void initialize(Config *config, Gas_<T> *gas) = 0;
    virtual void computeFlux(
      T const *leftVariable, T const *rightVariable, 
      double const *nonDimNormal, T *flux) = 0;
    
};

#include "../src/Flux.tpp"
typedef Flux_<double> Flux; // Instantiate with double data type for normal use.

#endif // FLUX_HPP 