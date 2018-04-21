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
#ifndef FLUX_SLIP_WALL_HPP
#define FLUX_SLIP_WALL_HPP

#include "../include/Flux.hpp"
#include "../include/GasPerfect.hpp"

template <typename T>
class FluxSlipWall_: public Flux_<T> {
    
  public:
  
    FluxSlipWall_();
    ~FluxSlipWall_();
  
    void initialize(Config *config, Gas_<T> *gas);
    void computeFlux(
      T const *leftVariable, T const *rightVariable, 
      double const *nonDimNormal, T *flux);
    
};

#include "../src/FluxSlipWall.tpp"
typedef FluxSlipWall_<double> FluxSlipWall;

#endif // FLUX_SLIP_WALL_HPP