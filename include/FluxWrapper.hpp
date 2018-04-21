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
#ifndef FLUX_WRAPPER_HPP
#define FLUX_WRAPPER_HPP

#include <Eigen/Dense>

#include "../include/Config.hpp"
#include "../include/GasPerfect.hpp"

template < 
  template<typename> typename F, short N_LEFT, short N_RIGHT, short N_OUT 
>
class FluxWrapper_ {
  
  protected:
    
    double const *nonDimNormal_;
    Config *config_;
    
  public:
    
    enum {
      LeftInputAtCompileTime = N_LEFT,
      RightInputAtCompileTime = N_RIGHT,
      InputsAtCompileTime = N_LEFT+N_RIGHT,
      ValuesAtCompileTime = N_OUT
    };
    typedef Eigen::Matrix<double, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<double, ValuesAtCompileTime, 1> ValueType;
  
    FluxWrapper_();
    ~FluxWrapper_();
    
    void initialize(Config *config);
    void setNonDimNormal(double const * nonDimNormal);
    
    template <typename E>
    void operator() (
      const Eigen::Matrix<E, InputsAtCompileTime, 1> & vIn, 
      Eigen::Matrix<E, ValuesAtCompileTime, 1> * vOut) const;
    
};

#include "../src/FluxWrapper.tpp"

#endif // FLUX_WRAPPER_HPP