#ifndef READ_IGLOO_FILES_H
#define READ_IGLOO_FILES_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

double NEO_integral_flux(double m);

void H_readInt_FromFile(ifstream& file, int& firstInt);

void H_getRowCol_FromFile(string fileName);

void readIgloo(string lonDirectory, string dens, int lat, vector<double>& iglooData);

void readVelDist(string vel_fn, vector<double>& velDist);

void setupNEOIgloo(string lonDirectory, int lat, vector<double>& velDist, vector<double>& iglooDataHi, vector<double>& iglooDataNEO, double m_min, double m_max);


#endif 