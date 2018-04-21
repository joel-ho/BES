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

#include "../include/FluxJacobianWrapper.hpp"

template <typename F>
FluxJacobianWrapper_<F>::FluxJacobianWrapper_ () {}

template <typename F>
FluxJacobianWrapper_<F>::~FluxJacobianWrapper_ () {}

template <typename F>
void FluxJacobianWrapper_<F>::initialize(Config *config) {
  config_ = config;
}

template <typename F>
void FluxJacobianWrapper_<F>::computeJacobian
(
  double const * leftVariable,
  double const * rightVariable,
  double const * nonDimNormal,
  double * leftJacobian,
  double * rightJacobian
) 
{
  
  // Required Eigen data structures
  double const dummyRightVariable[1]={0};
  double const *internalRightVariable;
  Eigen::Matrix<double, F::InputsAtCompileTime, 1> eigenInputVector;
  Eigen::Matrix<double, F::ValuesAtCompileTime, 1> eigenOutputVector;
  Eigen::Matrix<double, F::ValuesAtCompileTime, F::InputsAtCompileTime> eigenJacobian;
  
  // Template class autodiff function
  Eigen::AutoDiffJacobian<F> fAuto;
  
  // Map input from raw pointer array to eigen matrix
  internalRightVariable = (rightVariable != NULL)? rightVariable:dummyRightVariable;
  eigenInputVector << 
    Eigen::Map< const Eigen::Matrix<double, F::LeftInputAtCompileTime, 1> > 
      ( leftVariable, F::LeftInputAtCompileTime, 1 ), 
    Eigen::Map< const Eigen::Matrix<double, F::RightInputAtCompileTime, 1> > 
      ( internalRightVariable, F::RightInputAtCompileTime, 1 );
  
  // Solve for jacobian
  fAuto.initialize(config_);
  fAuto.setNonDimNormal(nonDimNormal);
  fAuto(eigenInputVector, &eigenOutputVector, &eigenJacobian);
  
  // Copy to raw pointer array
  for (short iOut=0; iOut<F::ValuesAtCompileTime; iOut++) {
    for (short iIn=0; iIn<F::LeftInputAtCompileTime; iIn++) {
      leftJacobian[ iOut*F::LeftInputAtCompileTime + iIn ] = 
        eigenJacobian(iOut, iIn);
    }
  }
  if (rightJacobian != NULL) {
    for (short iOut=0; iOut<F::ValuesAtCompileTime; iOut++) {
      for (short iIn=0; iIn<F::RightInputAtCompileTime; iIn++) {
        rightJacobian[ iOut*F::RightInputAtCompileTime + iIn ] =  
          eigenJacobian(iOut, iIn+F::LeftInputAtCompileTime);
      }
    }
  }
  
}
