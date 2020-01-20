#ifndef GRID_H
#define GRID_H

#include <iostream>

#include "Alloc.h"
#include "Parameters.h"
#include "PrecisionTypes.h"

/** Grid Data */
struct grid {
  /** number of cells - X direction, including + 2 (guard cells) */
  int nxc;
  /** number of nodes - X direction, including + 2 extra nodes for guard cells
   */
  int nxn;
  /** number of cell - Y direction, including + 2 (guard cells) */
  int nyc;
  /** number of nodes - Y direction, including + 2 extra nodes for guard cells
   */
  int nyn;
  /** number of cell - Z direction, including + 2 (guard cells) */
  int nzc;
  /** number of nodes - Z direction, including + 2 extra nodes for guard cells
   */
  int nzn;
  /** dx = space step - X direction */
  double dx;
  /** dy = space step - Y direction */
  double dy;
  /** dz = space step - Z direction */
  double dz;
  /** invdx = 1/dx */
  FPfield invdx;
  /** invdy = 1/dy */
  FPfield invdy;
  /** invdz = 1/dz */
  FPfield invdz;
  /** invol = inverse of volume*/
  FPfield invVOL;
  /** local grid boundaries coordinate  */
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;
  /** domain size */
  double Lx, Ly, Lz;

  /** Periodicity for fields X **/
  bool PERIODICX;
  /** Periodicity for fields Y **/
  bool PERIODICY;
  /** Periodicity for fields Z **/
  bool PERIODICZ;

  // Nodes coordinate
  /** coordinate node X */
  FPfield *XN_flat;
  FPfield ***XN;
  /** coordinate node Y */
  FPfield *YN_flat;
  FPfield ***YN;
  /** coordinate node Z */
  FPfield *ZN_flat;
  FPfield ***ZN;
};

/** Set up the grid quantities */
void setGrid(struct parameters *param, struct grid *grd);

/** Set up the grid quantities */
void printGrid(struct grid *grd);

/** allocate electric and magnetic field */
void grid_deallocate(struct grid *grd);

/** interpolation Node to Center */
void interpN2Cfield(FPfield ***vecFieldCx, FPfield ***vecFieldCy,
                    FPfield ***vecFieldCz, FPfield ***vecFieldNx,
                    FPfield ***vecFieldNy, FPfield ***vecFieldNz,
                    struct grid *grd);

/** interpolation Node to Center */
void interpC2Ninterp(FPinterp ***vecFieldN, FPinterp ***vecFieldC,
                     struct grid *grd);

/** interpolation Node to Center */
void interpC2Nfield(FPfield ***vecFieldN, FPfield ***vecFieldC,
                    struct grid *grd);

/** interpolation Node to Center */
void interpN2Cinterp(FPinterp ***vecFieldC, FPinterp ***vecFieldN,
                     struct grid *grd);

/** calculate gradient on nodes, given a scalar field defined on central points
 */
void gradC2N(FPfield ***gradXN, FPfield ***gradYN, FPfield ***gradZN,
             FPfield ***scFieldC, grid *grd);

/** calculate gradient on nodes, given a scalar field defined on central points
 */
void gradN2C(FPfield ***gradXC, FPfield ***gradYC, FPfield ***gradZC,
             FPfield ***scFieldN, grid *grd);

/** calculate divergence on central points, given a vector field defined on
 * nodes  */
void divN2C(FPfield ***divC, FPfield ***vecFieldXN, FPfield ***vecFieldYN,
            FPfield ***vecFieldZN, grid *grd);

/** calculate divergence on central points, given a Tensor field defined on
 * nodes  */
void divSymmTensorN2C(FPinterp ***divCX, FPinterp ***divCY, FPinterp ***divCZ,
                      FPinterp ***pXX, FPinterp ***pXY, FPinterp ***pXZ,
                      FPinterp ***pYY, FPinterp ***pYZ, FPinterp ***pZZ,
                      grid *grd);

/** calculate divergence on nodes, given a vector field defined on central
 * points  */
void divC2N(FPfield ***divN, FPfield ***vecFieldXC, FPfield ***vecFieldYC,
            FPfield ***vecFieldZC, grid *grd);

/** calculate curl on nodes, given a vector field defined on central points  */
void curlC2N(FPfield ***curlXN, FPfield ***curlYN, FPfield ***curlZN,
             FPfield ***vecFieldXC, FPfield ***vecFieldYC,
             FPfield ***vecFieldZC, grid *grd);

/** calculate curl on central points, given a vector field defined on nodes  */
void curlN2C(FPfield ***curlXC, FPfield ***curlYC, FPfield ***curlZC,
             FPfield ***vecFieldXN, FPfield ***vecFieldYN,
             FPfield ***vecFieldZN, grid *grd);

/** calculate laplacian on nodes, given a scalar field defined on nodes */
void lapN2N(FPfield ***lapN, FPfield ***scFieldN, grid *grd);

/** calculate laplacian on nodes, given a scalar field defined on nodes */
void lapN2N_V(FPfield ***lapN, FPfield ***scFieldN, grid *grd);

/** calculate laplacian on central points, given a scalar field defined on
 * central points */
void lapC2C(FPfield ***lapC, FPfield ***scFieldC, grid *grd);

#endif
