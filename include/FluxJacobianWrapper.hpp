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
#ifndef FLUX_JACOBIAN_WRAPPER_HPP
#define FLUX_JACOBIAN_WRAPPER_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>

#include "../include/Options.hpp"
#include "../include/Config.hpp"
#include "../include/FluxRoe.hpp"
#include "../include/FluxFreestream.hpp"
#include "../include/FluxSlipWall.hpp"
#include "../include/FluxWrapper.hpp"
#include "../include/FluxJacobian.hpp"

// Templated with FluxWrapper_ class F
template <typename F> 
class FluxJacobianWrapper_: public FluxJacobian {
  
  protected:
    
    Config *config_;
    
  public:
    
    FluxJacobianWrapper_();
    ~FluxJacobianWrapper_();
    
    void initialize(Config *config);
    void computeJacobian(
      double const * leftVariable,
      double const * rightVariable,
      double const * nonDimNormal,
      double * leftJacobian,
      double * rightJacobian);
    
};

#include "../src/FluxJacobianWrapper.tpp"

typedef FluxJacobianWrapper_<FluxWrapper_<FluxRoe_, 
  static_cast<short>(NumVariable::EULER_TWO_DIM_CONSERVATIVE), 
  static_cast<short>(NumVariable::EULER_TWO_DIM_CONSERVATIVE), 
  static_cast<short>(NumVariable::EULER_TWO_DIM_CONSERVATIVE) > 
  > FluxJacobianRoeTwoDim;
  
typedef FluxJacobianWrapper_<FluxWrapper_<FluxRoe_, 
  static_cast<short>(NumVariable::EULER_THREE_DIM_CONSERVATIVE), 
  static_cast<short>(NumVariable::EULER_THREE_DIM_CONSERVATIVE), 
  static_cast<short>(NumVariable::EULER_THREE_DIM_CONSERVATIVE) > 
  > FluxJacobianRoeThreeDim;
  
typedef FluxJacobianWrapper_<FluxWrapper_<FluxFreestream_, 
  static_cast<short>(NumVariable::EULER_TWO_DIM_CONSERVATIVE), 
  static_cast<short>(NumVariable::EULER_TWO_DIM_PRIMITIVE), 
  static_cast<short>(NumVariable::EULER_TWO_DIM_CONSERVATIVE) > 
  > FluxJacobianFreestreamTwoDim;
  
typedef FluxJacobianWrapper_<FluxWrapper_<FluxFreestream_, 
  static_cast<short>(NumVariable::EULER_THREE_DIM_CONSERVATIVE), 
  static_cast<short>(NumVariable::EULER_THREE_DIM_PRIMITIVE), 
  static_cast<short>(NumVariable::EULER_THREE_DIM_CONSERVATIVE) > 
  > FluxJacobianFreestreamThreeDim;
  
typedef FluxJacobianWrapper_<FluxWrapper_<FluxSlipWall_, 
  static_cast<short>(NumVariable::EULER_TWO_DIM_CONSERVATIVE), 
  1, 
  static_cast<short>(NumVariable::EULER_TWO_DIM_CONSERVATIVE) > 
  > FluxJacobianSlipWallTwoDim;
  
typedef FluxJacobianWrapper_<FluxWrapper_<FluxSlipWall_, 
  static_cast<short>(NumVariable::EULER_THREE_DIM_CONSERVATIVE), 
  1, 
  static_cast<short>(NumVariable::EULER_THREE_DIM_CONSERVATIVE) > 
  > FluxJacobianSlipWallThreeDim;

#endif // FLUX_JACOBIAN_WRAPPER_HPP