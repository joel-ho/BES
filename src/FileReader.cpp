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

#include "../include/FileReader.hpp"

void ReadFile::splitString (
  const string &stringIn, const string &delimiter, string *vectorOut) {

  int i=0, start=0, end;
  string trimmedStr;
  
  end = stringIn.find(delimiter);
  while (end != string::npos) {
    trimmedStr = stringIn.substr(start, end - start);
    trimWhiteSpace(trimmedStr);
    vectorOut[i] = trimmedStr;
    start = end + delimiter.length();
    end = stringIn.find(delimiter, start);
    i++;
  }
  trimmedStr = stringIn.substr(start, end);
  trimWhiteSpace(trimmedStr);
  vectorOut[i] = trimmedStr;
  
}

void ReadFile::trimWhiteSpace(string &stringIn) {
  
  string::iterator endPos;
  
  endPos = remove(stringIn.begin(), stringIn.end(), ' ');
  stringIn.erase(endPos, stringIn.end());
  
  endPos = remove(stringIn.begin(), stringIn.end(), '\t');
  stringIn.erase(endPos, stringIn.end());
  
  endPos = remove(stringIn.begin(), stringIn.end(), '\r');
  stringIn.erase(endPos, stringIn.end());
  
  endPos = remove(stringIn.begin(), stringIn.end(), '\0');
  stringIn.erase(endPos, stringIn.end());
  
}

void ReadFile::splitStringToNum (
  const string &stringIn, const string &delimiter, double *vectorOut) {
  
  int i=0, start=0, end;
  string trimmedStr;
  
  end = stringIn.find(delimiter);
  while (end != string::npos) {
    trimmedStr = stringIn.substr(start, end - start);
    trimWhiteSpace(trimmedStr);
    vectorOut[i] = stod(trimmedStr);
    start = end + delimiter.length();
    end = stringIn.find(delimiter, start);
    i++;
  }
  trimmedStr = stringIn.substr(start, end);
  trimWhiteSpace(trimmedStr);
  vectorOut[i] = stod(trimmedStr);
  
}

void ReadFile::encodeFaceIdentifier (
  string &identifier, vector<unsigned long> vertexNumbers) {
  // Create unique identifier for face given face vertex numbers
  
  stringstream identifierStream;
  
  sort(vertexNumbers.begin(), vertexNumbers.end());
  copy(vertexNumbers.begin(), vertexNumbers.end(),
    ostream_iterator<long>(identifierStream, "-"));
  identifier = identifierStream.str();
  identifier = identifier.substr(0, identifier.length()-1);
  
}
