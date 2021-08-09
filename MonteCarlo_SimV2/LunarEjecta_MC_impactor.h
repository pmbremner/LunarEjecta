// Impactor information
#ifndef LUNAREJECTA_MC_IMPACTOR_H
#define LUNAREJECTA_MC_IMPACTOR_H


#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <time.h>
#include <algorithm> // max(a,b)
//#include <mpi.h>
#include <stdlib.h> // rand(), atof()

using namespace std;



// Read igloo file from disk
//
// Input: igloo file name
//
// Output: igloo data
//         speed array (km/s)
//         I   J   PHI1   PHI2 THETA1 THETA2 PHIavg THETAavg


// Compute CDF from igloo file
//
// Input: igloo data
//
// Output: CDF
//         Normalization factor (#/m^2/yr)


// Choose an impactor randomly
//
// Input: Location on Moon
//
// Output: type of impactor (MEM_lo, MEM_hi, NEO)
//         speed (km/s)
//         horizon angle (degrees > 0)
//         azimuth angle (degrees, 0 = East, 90 = North)
//
// Need: CDF of igloo, and array of speed, horz, and azm ang




#endif 