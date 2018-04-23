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

#define BES_VERSION 0.1

#include <iostream>
#include <fstream>
#include <string>
#include <exception>

#include "../include/Options.hpp"
#include "../include/Config.hpp"
#include "../include/GasPerfect.hpp"
#include "../include/Mesh.hpp"
#include "../include/Solution.hpp"
#include "../include/GradientGreenGauss.hpp"
#include "../include/LimiterVenkatakrishnan.hpp"
#include "../include/FluxRoe.hpp"
#include "../include/FluxFreestream.hpp"
#include "../include/SolverExplicit.hpp"
#include "../include/SolverImplicit.hpp"

using namespace std;

int main (int argc, char *argv[]) {
  
  string configFilePath;
  
  Config *config;
  Gas *gas;
  Mesh *mesh;
  Solver *solver;
  
  cout << endl;
  cout << "###################################################################" << endl;
  cout << "###################################################################" << endl;
  cout << "" << endl;
  cout << "   ###############      ###################         ###########   " << endl;
  cout << "   ################     ###################       #############" << endl;
  cout << "      ###        ###      ###                    ###" << endl;
  cout << "      ###          ###    ###                  ###   " << endl;
  cout << "      ###           ###   ###                 ###" << endl;
  cout << "      ###           ###   ###                 ###" << endl;
  cout << "      ###          ###    ###                  ###" << endl;
  cout << "      ###        ####     ###                   ###" << endl;
  cout << "      #############       ################        #######" << endl;
  cout << "      #############       ################          ####### " << endl;
  cout << "      ###        ###      ###                             ###" << endl;
  cout << "      ###          ###    ###                               ###" << endl;
  cout << "      ###           ###   ###                                ###" << endl;
  cout << "      ###           ###   ###                                ###      " << endl;
  cout << "      ###           ###   ###                                ###" << endl;
  cout << "      ###          ###    ###                               ###" << endl;
  cout << "      ###        ###      ###                             ###" << endl;
  cout << "   ################      ##################    #############" << endl;
  cout << "   ###############       ##################    ###########" << endl;
  cout << "" << endl;
  cout << "   Basic Euler Solver v" << BES_VERSION << endl;
  cout << "   Created by: Ho Mun Onn, Joel" << endl;
  cout << "###################################################################" << endl;
  cout << "###################################################################" << endl << endl;
  
  if (argc != 2) {
    cout << "ERROR: Incorrect number of arguments given. Exiting..." << endl;
    exit(EXIT_FAILURE);
  }
  configFilePath = argv[1];
  cout << "BES input file: \"" << configFilePath << "\"." << endl << endl;
  
  config = new Config();
	mesh = new Mesh();
  gas = new GasPerfect();
  
  config->initialize(configFilePath);
  mesh->initialize(config);
  config->initializePostMesh(mesh);
  gas->initialize(config);
  
  switch (config->getSolverType()) {
  
    case SolverType::EXPLICIT:
      solver = new SolverExplicit();
      break;
      
    case SolverType::IMPLICIT:
      solver = new SolverImplicit();
      break;
      
    default:
      break;
  
  }
  
  solver->initialize(config, gas, mesh);
  solver->setUniversalInitialCondition(config->getFreestreamPrimitiveVariables(0));
  solver->solve();
  
  cout << endl << "Cleaning up... ";
  delete solver;
  delete gas;
  delete mesh;
  delete config;
  cout << "Done." << endl;
  
  return EXIT_SUCCESS;
}
