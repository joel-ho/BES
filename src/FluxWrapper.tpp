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

#include "../include/FluxWrapper.hpp"

template < template<typename> typename F, short N_LEFT, short N_RIGHT, short N_OUT >
FluxWrapper_<F, N_LEFT, N_RIGHT, N_OUT>::FluxWrapper_ () {}

template < template<typename> typename F, short N_LEFT, short N_RIGHT, short N_OUT >
FluxWrapper_<F, N_LEFT, N_RIGHT, N_OUT>::~FluxWrapper_ () {}

template < template<typename> typename F, short N_LEFT, short N_RIGHT, short N_OUT >
void FluxWrapper_<F, N_LEFT, N_RIGHT, N_OUT>::initialize (Config *config) {
  config_ = config;
}

template < template<typename> typename F, short N_LEFT, short N_RIGHT, short N_OUT >
void FluxWrapper_<F, N_LEFT, N_RIGHT, N_OUT>::setNonDimNormal(double const * nonDimNormal) {
  nonDimNormal_ = nonDimNormal;
}

template < template<typename> typename F, short N_LEFT, short N_RIGHT, short N_OUT >
template<typename E>
void FluxWrapper_<F, N_LEFT, N_RIGHT, N_OUT>::operator() (
  Eigen::Matrix<E, InputsAtCompileTime, 1> const & vIn, 
  Eigen::Matrix<E, ValuesAtCompileTime, 1> * vOut
) const
{
  
  GasPerfect_<E> gas;
  F<E> f;
  
  gas.initialize(config_);
  f.initialize(config_, &gas);
  
  f.computeFlux(
    &(vIn[0]),
    &(vIn[N_LEFT]),
    nonDimNormal_,
    &((*vOut)[0])
  );
  
}