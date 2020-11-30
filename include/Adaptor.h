#ifndef FEADAPTOR_HEADER
#define FEADAPTOR_HEADER

// For iPic3D arrays
#include "Alloc.h"
// Access to physical quantities
#include "EMfield.h"

namespace Adaptor {
void Initialize(const int ns, const double B0x, const double B0y, const double B0z,
                const int start_x, const int start_y, const int start_z,
                const int nx, const int ny, const int nz,
                const double dx, const double dy, const double dz,
                const char *paraview_script_path);

void Finalize();

void CoProcess(double time, unsigned int timeStep, EMfield *EMf);
} // namespace Adaptor

#endif
