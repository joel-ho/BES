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
#ifndef FLUX_JACOBIAN_HPP
#define FLUX_JACOBIAN_HPP

#include "../include/Config.hpp"
#include "../include/CoefficientMatrix.hpp"

class FluxJacobian {
  
  protected:
    Config *config_;
    
  public:
    
    FluxJacobian();
    virtual ~FluxJacobian();
    
    virtual void initialize(Config *config) = 0;
    
    // Jacobians given in row major order
    virtual void computeJacobian(
      double const * leftVariable,
      double const * rightVariable,
      double const * nonDimNormal,
      double * leftJacobian,
      double * rightJacobian) = 0;
    
};

#endif // FLUX_JACOBIAN_HPP