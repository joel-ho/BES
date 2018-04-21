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
#ifndef LIMITER_VENKATAKRISHNAN_HPP
#define LIMITER_VENKATAKRISHNAN_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <vector>
#include <algorithm>

#include "../include/Limiter.hpp"
#include "../include/GradientGreenGauss.hpp"

class LimiterVenkatakrishnan: public Limiter {
  private:
    double computeFaceLimiter_(double &deltaOneMin, double &deltaOneMax, 
      double &deltaTwo, double &epsilon);
  public:
    void computeLimiterTerm(
      SolutionSnapshot *solutionSnapshot, Gradient *gradient);
      
};

#endif // LIMITER_VENKATAKRISHNAN_HPP