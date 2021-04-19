#ifndef ADAPTOR_H
#define ADAPTOR_H

// For iPic3D arrays
#include "PrecisionTypes.h"

namespace Adaptor
{
void Initialize(const char* script, const int start_x, const int start_y, const int start_z, \
                          const int nx, const int ny, const int nz, \
                          const double dx, const double dy, const double dz);

void Finalize();

void CoProcess(
  double time, unsigned int timeStep, FPfield ***Bx, FPfield ***By, FPfield ***Bz, FPinterp ***rhons0_array, FPinterp ***rhons1_array);
}

#endif
