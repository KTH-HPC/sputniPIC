#include "BC.h"

/** Put Boundary Conditions on Boundaries */

//////////
// POPULATE GHOST CELL ON NODES
//////////

/** Apply BC to scalar interp quantity defined on nodes - Interpolation quantity
 */
void applyBCscalarDensN(FPinterp ***scalarN, grid *grd, parameters *param) {
  ///////////////////////
  ///
  ///    FACE
  ///

  // X direction
  if (param->PERIODICX == true) {
    for (int j = 1; j < grd->nyn - 1; j++)
      for (int k = 1; k < grd->nzn - 1; k++) {
        scalarN[0][j][k] = scalarN[grd->nxn - 3][j][k];
        scalarN[grd->nxn - 1][j][k] = scalarN[2][j][k];
      }
  }

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == true) {
    for (int i = 1; i < grd->nxn - 1; i++)
      for (int k = 1; k < grd->nzn - 1; k++) {
        // rhon
        scalarN[i][0][k] = scalarN[i][grd->nyn - 3][k];
        scalarN[i][grd->nyn - 1][k] = scalarN[i][2][k];
      }
  }

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == true) {
    for (int i = 1; i < grd->nxn - 1; i++)
      for (int j = 1; j < grd->nyn - 1; j++) {
        // rhon
        scalarN[i][j][0] = scalarN[i][j][grd->nzn - 3];
        scalarN[i][j][grd->nzn - 1] = scalarN[i][j][2];
      }
  }

  ///////////////////////
  ///
  ///    EDGES
  ///

  // X-EDGE
  if (param->PERIODICY == true || param->PERIODICZ == true) {
    for (int i = 1; i < (grd->nxn - 1); i++) {
      scalarN[i][grd->nyn - 1][grd->nzn - 1] = scalarN[i][2][2];
      scalarN[i][0][0] = scalarN[i][grd->nyn - 3][grd->nzn - 3];
      scalarN[i][0][grd->nzn - 1] = scalarN[i][grd->nyn - 3][2];
      scalarN[i][grd->nyn - 1][0] = scalarN[i][2][grd->nzn - 3];
    }
  }

  // Y-EDGE
  if (param->PERIODICX == true || param->PERIODICZ == true) {
    for (int i = 1; i < (grd->nyn - 1); i++) {
      scalarN[grd->nxn - 1][i][grd->nzn - 1] = scalarN[2][i][2];
      scalarN[0][i][0] = scalarN[grd->nxn - 3][i][grd->nzn - 3];
      scalarN[0][i][grd->nzn - 1] = scalarN[grd->nxn - 3][i][2];
      scalarN[grd->nxn - 1][i][0] = scalarN[2][i][grd->nzn - 3];
    }
  }

  // Z-EDGE
  if (param->PERIODICX == true || param->PERIODICY == true) {
    for (int i = 1; i < (grd->nzn - 1); i++) {
      scalarN[grd->nxn - 1][grd->nyn - 1][i] = scalarN[2][2][i];
      scalarN[0][0][i] = scalarN[grd->nxn - 3][grd->nyn - 3][i];
      scalarN[grd->nxn - 1][0][i] = scalarN[2][grd->nyn - 3][i];
      scalarN[0][grd->nyn - 1][i] = scalarN[grd->nxn - 3][2][i];
    }
  }

  // Corners
  if (param->PERIODICX == true || param->PERIODICY == true ||
      param->PERIODICZ == true) {
    scalarN[grd->nxn - 1][grd->nyn - 1][grd->nzn - 1] = scalarN[2][2][2];
    scalarN[0][grd->nyn - 1][grd->nzn - 1] = scalarN[grd->nxn - 3][2][2];
    scalarN[grd->nxn - 1][0][grd->nzn - 1] = scalarN[2][grd->nyn - 3][2];
    scalarN[0][0][grd->nzn - 1] = scalarN[grd->nxn - 3][grd->nyn - 3][2];
    scalarN[grd->nxn - 1][grd->nyn - 1][0] = scalarN[2][2][grd->nzn - 3];
    scalarN[0][grd->nyn - 1][0] = scalarN[grd->nxn - 3][2][grd->nzn - 3];
    scalarN[grd->nxn - 1][0][0] = scalarN[2][grd->nyn - 3][grd->nzn - 3];
    scalarN[0][0][0] = scalarN[grd->nxn - 3][grd->nyn - 3][grd->nzn - 3];
  }

  // FACE NON PERIODIC
  // X direction
  if (param->PERIODICX == false) {
    for (int j = 0; j < grd->nyn; j++)
      for (int k = 0; k < grd->nzn; k++) {
        scalarN[0][j][k] = scalarN[1][j][k];
        scalarN[grd->nxn - 1][j][k] = scalarN[grd->nxn - 2][j][k];
      }
  }

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == false) {
    for (int i = 0; i < grd->nxn; i++)
      for (int k = 0; k < grd->nzn; k++) {
        scalarN[i][0][k] = scalarN[i][1][k];
        scalarN[i][grd->nyn - 1][k] = scalarN[i][grd->nyn - 2][k];
      }
  }

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == false) {
    for (int i = 0; i < grd->nxn; i++)
      for (int j = 0; j < grd->nyn; j++) {
        scalarN[i][j][0] = scalarN[i][j][1];
        scalarN[i][j][grd->nzn - 1] = scalarN[i][j][grd->nzn - 2];
      }
  }
}

/** Apply BC to scalar interp quantity defined on nodes - Interpolation quantity
 */
void applyBCscalarFieldN(FPfield ***scalarN, grid *grd, parameters *param) {
  ///////////////////////
  ///
  ///    FACE
  ///

  // X direction
  if (param->PERIODICX == true) {
    for (int j = 1; j < grd->nyn - 1; j++)
      for (int k = 1; k < grd->nzn - 1; k++) {
        scalarN[0][j][k] = scalarN[grd->nxn - 3][j][k];
        scalarN[grd->nxn - 1][j][k] = scalarN[2][j][k];
      }
  }

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == true) {
    for (int i = 1; i < grd->nxn - 1; i++)
      for (int k = 1; k < grd->nzn - 1; k++) {
        // rhon
        scalarN[i][0][k] = scalarN[i][grd->nyn - 3][k];
        scalarN[i][grd->nyn - 1][k] = scalarN[i][2][k];
      }
  }

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == true) {
    for (int i = 1; i < grd->nxn - 1; i++)
      for (int j = 1; j < grd->nyn - 1; j++) {
        // rhon
        scalarN[i][j][0] = scalarN[i][j][grd->nzn - 3];
        scalarN[i][j][grd->nzn - 1] = scalarN[i][j][2];
      }
  }

  ///////////////////////
  ///
  ///    EDGES
  ///

  // X-EDGE
  if (param->PERIODICY == true || param->PERIODICZ == true) {
    for (int i = 1; i < (grd->nxn - 1); i++) {
      scalarN[i][grd->nyn - 1][grd->nzn - 1] = scalarN[i][2][2];
      scalarN[i][0][0] = scalarN[i][grd->nyn - 3][grd->nzn - 3];
      scalarN[i][0][grd->nzn - 1] = scalarN[i][grd->nyn - 3][2];
      scalarN[i][grd->nyn - 1][0] = scalarN[i][2][grd->nzn - 3];
    }
  }

  // Y-EDGE
  if (param->PERIODICX == true || param->PERIODICZ == true) {
    for (int i = 1; i < (grd->nyn - 1); i++) {
      scalarN[grd->nxn - 1][i][grd->nzn - 1] = scalarN[2][i][2];
      scalarN[0][i][0] = scalarN[grd->nxn - 3][i][grd->nzn - 3];
      scalarN[0][i][grd->nzn - 1] = scalarN[grd->nxn - 3][i][2];
      scalarN[grd->nxn - 1][i][0] = scalarN[2][i][grd->nzn - 3];
    }
  }

  // Z-EDGE
  if (param->PERIODICX == true || param->PERIODICY == true) {
    for (int i = 1; i < (grd->nzn - 1); i++) {
      scalarN[grd->nxn - 1][grd->nyn - 1][i] = scalarN[2][2][i];
      scalarN[0][0][i] = scalarN[grd->nxn - 3][grd->nyn - 3][i];
      scalarN[grd->nxn - 1][0][i] = scalarN[2][grd->nyn - 3][i];
      scalarN[0][grd->nyn - 1][i] = scalarN[grd->nxn - 3][2][i];
    }
  }

  // Corners
  if (param->PERIODICX == true || param->PERIODICY == true ||
      param->PERIODICZ == true) {
    scalarN[grd->nxn - 1][grd->nyn - 1][grd->nzn - 1] = scalarN[2][2][2];
    scalarN[0][grd->nyn - 1][grd->nzn - 1] = scalarN[grd->nxn - 3][2][2];
    scalarN[grd->nxn - 1][0][grd->nzn - 1] = scalarN[2][grd->nyn - 3][2];
    scalarN[0][0][grd->nzn - 1] = scalarN[grd->nxn - 3][grd->nyn - 3][2];
    scalarN[grd->nxn - 1][grd->nyn - 1][0] = scalarN[2][2][grd->nzn - 3];
    scalarN[0][grd->nyn - 1][0] = scalarN[grd->nxn - 3][2][grd->nzn - 3];
    scalarN[grd->nxn - 1][0][0] = scalarN[2][grd->nyn - 3][grd->nzn - 3];
    scalarN[0][0][0] = scalarN[grd->nxn - 3][grd->nyn - 3][grd->nzn - 3];
  }

  // FACE NON PERIODIC
  // X direction
  if (param->PERIODICX == false) {
    for (int j = 0; j < grd->nyn; j++)
      for (int k = 0; k < grd->nzn; k++) {
        scalarN[0][j][k] = scalarN[1][j][k];
        scalarN[grd->nxn - 1][j][k] = scalarN[grd->nxn - 2][j][k];
      }
  }

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == false) {
    for (int i = 0; i < grd->nxn; i++)
      for (int k = 0; k < grd->nzn; k++) {
        scalarN[i][0][k] = scalarN[i][1][k];
        scalarN[i][grd->nyn - 1][k] = scalarN[i][grd->nyn - 2][k];
      }
  }

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == false) {
    for (int i = 0; i < grd->nxn; i++)
      for (int j = 0; j < grd->nyn; j++) {
        scalarN[i][j][0] = scalarN[i][j][1];
        scalarN[i][j][grd->nzn - 1] = scalarN[i][j][grd->nzn - 2];
      }
  }
}

///////// USE THIS TO IMPOSE BC TO ELECTRIC FIELD
///////// NOW THIS IS FIXED TO ZERO

/** Apply BC to scalar interp quantity defined on nodes - Interpolation quantity
 */
void applyBCscalarFieldNzero(FPfield ***scalarN, grid *grd, parameters *param) {
  ///////////////////////
  ///
  ///    FACE
  ///

  // X direction
  if (param->PERIODICX == true) {
    for (int j = 1; j < grd->nyn - 1; j++)
      for (int k = 1; k < grd->nzn - 1; k++) {
        scalarN[0][j][k] = scalarN[grd->nxn - 3][j][k];
        scalarN[grd->nxn - 1][j][k] = scalarN[2][j][k];
      }
  }

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == true) {
    for (int i = 1; i < grd->nxn - 1; i++)
      for (int k = 1; k < grd->nzn - 1; k++) {
        // rhon
        scalarN[i][0][k] = scalarN[i][grd->nyn - 3][k];
        scalarN[i][grd->nyn - 1][k] = scalarN[i][2][k];
      }
  }

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == true) {
    for (int i = 1; i < grd->nxn - 1; i++)
      for (int j = 1; j < grd->nyn - 1; j++) {
        // rhon
        scalarN[i][j][0] = scalarN[i][j][grd->nzn - 3];
        scalarN[i][j][grd->nzn - 1] = scalarN[i][j][2];
      }
  }

  ///////////////////////
  ///
  ///    EDGES
  ///

  // X-EDGE
  if (param->PERIODICY == true || param->PERIODICZ == true) {
    for (int i = 1; i < (grd->nxn - 1); i++) {
      scalarN[i][grd->nyn - 1][grd->nzn - 1] = scalarN[i][2][2];
      scalarN[i][0][0] = scalarN[i][grd->nyn - 3][grd->nzn - 3];
      scalarN[i][0][grd->nzn - 1] = scalarN[i][grd->nyn - 3][2];
      scalarN[i][grd->nyn - 1][0] = scalarN[i][2][grd->nzn - 3];
    }
  }

  // Y-EDGE
  if (param->PERIODICX == true || param->PERIODICZ == true) {
    for (int i = 1; i < (grd->nyn - 1); i++) {
      scalarN[grd->nxn - 1][i][grd->nzn - 1] = scalarN[2][i][2];
      scalarN[0][i][0] = scalarN[grd->nxn - 3][i][grd->nzn - 3];
      scalarN[0][i][grd->nzn - 1] = scalarN[grd->nxn - 3][i][2];
      scalarN[grd->nxn - 1][i][0] = scalarN[2][i][grd->nzn - 3];
    }
  }

  // Z-EDGE
  if (param->PERIODICX == true || param->PERIODICY == true) {
    for (int i = 1; i < (grd->nzn - 1); i++) {
      scalarN[grd->nxn - 1][grd->nyn - 1][i] = scalarN[2][2][i];
      scalarN[0][0][i] = scalarN[grd->nxn - 3][grd->nyn - 3][i];
      scalarN[grd->nxn - 1][0][i] = scalarN[2][grd->nyn - 3][i];
      scalarN[0][grd->nyn - 1][i] = scalarN[grd->nxn - 3][2][i];
    }
  }

  // Corners
  if (param->PERIODICX == true || param->PERIODICY == true ||
      param->PERIODICZ == true) {
    scalarN[grd->nxn - 1][grd->nyn - 1][grd->nzn - 1] = scalarN[2][2][2];
    scalarN[0][grd->nyn - 1][grd->nzn - 1] = scalarN[grd->nxn - 3][2][2];
    scalarN[grd->nxn - 1][0][grd->nzn - 1] = scalarN[2][grd->nyn - 3][2];
    scalarN[0][0][grd->nzn - 1] = scalarN[grd->nxn - 3][grd->nyn - 3][2];
    scalarN[grd->nxn - 1][grd->nyn - 1][0] = scalarN[2][2][grd->nzn - 3];
    scalarN[0][grd->nyn - 1][0] = scalarN[grd->nxn - 3][2][grd->nzn - 3];
    scalarN[grd->nxn - 1][0][0] = scalarN[2][grd->nyn - 3][grd->nzn - 3];
    scalarN[0][0][0] = scalarN[grd->nxn - 3][grd->nyn - 3][grd->nzn - 3];
  }

  // FACE NON PERIODIC
  // X direction
  if (param->PERIODICX == false) {
    for (int j = 0; j < grd->nyn; j++)
      for (int k = 0; k < grd->nzn; k++) {
        scalarN[0][j][k] = 0.0;
        scalarN[grd->nxn - 1][j][k] = 0.0;
      }
  }

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == false) {
    for (int i = 0; i < grd->nxn; i++)
      for (int k = 0; k < grd->nzn; k++) {
        scalarN[i][0][k] = 0.0;
        scalarN[i][grd->nyn - 1][k] = 0.0;
      }
  }

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == false) {
    for (int i = 0; i < grd->nxn; i++)
      for (int j = 0; j < grd->nyn; j++) {
        scalarN[i][j][0] = 0.0;
        scalarN[i][j][grd->nzn - 1] = 0.0;
      }
  }
}

///////////////
////
////    add Densities
////
////
///////////////

// apply boundary conditions to species interpolated densities
void applyBCids(struct interpDensSpecies *ids, struct grid *grd,
                struct parameters *param) {
  /////////////////(///
  // apply BC on X
  /////

  // X direction
  if (param->PERIODICX == true) {
    for (int j = 1; j < grd->nyn - 1; j++)
      for (int k = 1; k < grd->nzn - 1; k++) {
        // rhon
        ids->rhon[1][j][k] += ids->rhon[grd->nxn - 2][j][k];
        ids->rhon[grd->nxn - 2][j][k] =
            ids->rhon[1][j][k];  // second is = not +=

        // Jx
        ids->Jx[1][j][k] += ids->Jx[grd->nxn - 2][j][k];
        ids->Jx[grd->nxn - 2][j][k] = ids->Jx[1][j][k];  // second is = not +=
        // Jy
        ids->Jy[1][j][k] += ids->Jy[grd->nxn - 2][j][k];
        ids->Jy[grd->nxn - 2][j][k] = ids->Jy[1][j][k];  // second is = not +=
        // Jz
        ids->Jz[1][j][k] += ids->Jz[grd->nxn - 2][j][k];
        ids->Jz[grd->nxn - 2][j][k] = ids->Jz[1][j][k];  // second is = not +=

        // pxx
        ids->pxx[1][j][k] += ids->pxx[grd->nxn - 2][j][k];
        ids->pxx[grd->nxn - 2][j][k] = ids->pxx[1][j][k];  // second is = not +=
        // pxy
        ids->pxy[1][j][k] += ids->pxy[grd->nxn - 2][j][k];
        ids->pxy[grd->nxn - 2][j][k] = ids->pxy[1][j][k];  // second is = not +=
        // pxz
        ids->pxz[1][j][k] += ids->pxz[grd->nxn - 2][j][k];
        ids->pxz[grd->nxn - 2][j][k] = ids->pxz[1][j][k];  // second is = not +=

        // pyy
        ids->pyy[1][j][k] += ids->pyy[grd->nxn - 2][j][k];
        ids->pyy[grd->nxn - 2][j][k] = ids->pyy[1][j][k];  // second is = not +=
        // pyz
        ids->pyz[1][j][k] += ids->pyz[grd->nxn - 2][j][k];
        ids->pyz[grd->nxn - 2][j][k] = ids->pyz[1][j][k];  // second is = not +=
        // pzz
        ids->pzz[1][j][k] += ids->pzz[grd->nxn - 2][j][k];
        ids->pzz[grd->nxn - 2][j][k] = ids->pzz[1][j][k];  // second is = not +=

      }  // end of loop over the grid
  }      // end of periodic in X direction

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == true) {
    for (int i = 1; i < grd->nxn - 1; i++)
      for (int k = 1; k < grd->nzn - 1; k++) {
        // rhon
        ids->rhon[i][1][k] += ids->rhon[i][grd->nyn - 2][k];
        ids->rhon[i][grd->nyn - 2][k] =
            ids->rhon[i][1][k];  // second is = not +=

        // Jx
        ids->Jx[i][1][k] += ids->Jx[i][grd->nyn - 2][k];
        ids->Jx[i][grd->nyn - 2][k] = ids->Jx[i][1][k];  // second is = not +=
        // Jy
        ids->Jy[i][1][k] += ids->Jy[i][grd->nyn - 2][k];
        ids->Jy[i][grd->nyn - 2][k] = ids->Jy[i][1][k];  // second is = not +=
        // Jz
        ids->Jz[i][1][k] += ids->Jz[i][grd->nyn - 2][k];
        ids->Jz[i][grd->nyn - 2][k] = ids->Jz[i][1][k];  // second is = not +=

        // pxx
        ids->pxx[i][1][k] += ids->pxx[i][grd->nyn - 2][k];
        ids->pxx[i][grd->nyn - 2][k] = ids->pxx[i][1][k];  // second is = not +=
        // pxy
        ids->pxy[i][1][k] += ids->pxy[i][grd->nyn - 2][k];
        ids->pxy[i][grd->nyn - 2][k] = ids->pxy[i][1][k];  // second is = not +=
        // pxz
        ids->pxz[i][1][k] += ids->pxz[i][grd->nyn - 2][k];
        ids->pxz[i][grd->nyn - 2][k] = ids->pxz[i][1][k];  // second is = not +=

        // pyy
        ids->pyy[i][1][k] += ids->pyy[i][grd->nyn - 2][k];
        ids->pyy[i][grd->nyn - 2][k] = ids->pyy[i][1][k];  // second is = not +=
        // pyz
        ids->pyz[i][1][k] += ids->pyz[i][grd->nyn - 2][k];
        ids->pyz[i][grd->nyn - 2][k] = ids->pyz[i][1][k];  // second is = not +=
        // pzz
        ids->pzz[i][1][k] += ids->pzz[i][grd->nyn - 2][k];
        ids->pzz[i][grd->nyn - 2][k] = ids->pzz[i][1][k];  // second is = not +=

      }  // end of loop

  }  // end of PERIODICY

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == true) {
    for (int i = 1; i < grd->nxn - 1; i++)
      for (int j = 1; j < grd->nyn - 1; j++) {
        // rhon
        ids->rhon[i][j][1] += ids->rhon[i][j][grd->nzn - 2];
        ids->rhon[i][j][grd->nzn - 2] =
            ids->rhon[i][j][1];  // second is = not +=

        // Jx
        ids->Jx[i][j][1] += ids->Jx[i][j][grd->nzn - 2];
        ids->Jx[i][j][grd->nzn - 2] = ids->Jx[i][j][1];  // second is = not +=
        // Jy
        ids->Jy[i][j][1] += ids->Jy[i][j][grd->nzn - 2];
        ids->Jy[i][j][grd->nzn - 2] = ids->Jy[i][j][1];  // second is = not +=
        // Jz
        ids->Jz[i][j][1] += ids->Jz[i][j][grd->nzn - 2];
        ids->Jz[i][j][grd->nzn - 2] = ids->Jz[i][j][1];  // second is = not +=

        // pxx
        ids->pxx[i][j][1] += ids->pxx[i][j][grd->nzn - 2];
        ids->pxx[i][j][grd->nzn - 2] = ids->pxx[i][j][1];  // second is = not +=
        // pxy
        ids->pxy[i][j][1] += ids->pxy[i][j][grd->nzn - 2];
        ids->pxy[i][j][grd->nzn - 2] = ids->pxy[i][j][1];  // second is = not +=
        // pxz
        ids->pxz[i][j][1] += ids->pxz[i][j][grd->nzn - 2];
        ids->pxz[i][j][grd->nzn - 2] = ids->pxz[i][j][1];  // second is = not +=
        // pyy
        ids->pyy[i][j][1] += ids->pyy[i][j][grd->nzn - 2];
        ids->pyy[i][j][grd->nzn - 2] = ids->pyy[i][j][1];  // second is = not +=
        // pyz
        ids->pyz[i][j][1] += ids->pyz[i][j][grd->nzn - 2];
        ids->pyz[i][j][grd->nzn - 2] = ids->pyz[i][j][1];  // second is = not +=
        // pzz
        ids->pzz[i][j][1] += ids->pzz[i][j][grd->nzn - 2];
        ids->pzz[i][j][grd->nzn - 2] = ids->pzz[i][j][1];  // second is = not +=
      }
  }

  // apply BC if BC are not periodic
  // X direction
  if (param->PERIODICX == false) {
    for (int j = 1; j < grd->nyn - 1; j++)
      for (int k = 1; k < grd->nzn - 1; k++) {
        // rhon
        ids->rhon[1][j][k] *= 2;
        ids->rhon[grd->nxn - 2][j][k] *= 2;  // second is = not +=

        // Jx
        ids->Jx[1][j][k] *= 2;
        ids->Jx[grd->nxn - 2][j][k] *= 2;
        // Jy
        ids->Jy[1][j][k] *= 2;
        ids->Jy[grd->nxn - 2][j][k] *= 2;
        // Jz
        ids->Jz[1][j][k] *= 2;
        ids->Jz[grd->nxn - 2][j][k] *= 2;

        // pxx
        ids->pxx[1][j][k] *= 2;
        ids->pxx[grd->nxn - 2][j][k] *= 2;
        // pxy
        ids->pxy[1][j][k] *= 2;
        ids->pxy[grd->nxn - 2][j][k] *= 2;
        // pxz
        ids->pxz[1][j][k] *= 2;
        ids->pxz[grd->nxn - 2][j][k] *= 2;

        // pyy
        ids->pyy[1][j][k] *= 2;
        ids->pyy[grd->nxn - 2][j][k] *= 2;
        // pyz
        ids->pyz[1][j][k] *= 2;
        ids->pyz[grd->nxn - 2][j][k] *= 2;
        // pzz
        ids->pzz[1][j][k] *= 2;
        ids->pzz[grd->nxn - 2][j][k] *= 2;

      }  // end of loop over the grid
  }      // end of not periodic in X direction

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == false) {
    for (int i = 1; i < grd->nxn - 1; i++)
      for (int k = 1; k < grd->nzn - 1; k++) {
        // rhon
        ids->rhon[i][1][k] *= 2;
        ids->rhon[i][grd->nyn - 2][k] *= 2;

        // Jx
        ids->Jx[i][1][k] *= 2;
        ids->Jx[i][grd->nyn - 2][k] *= 2;
        // Jy
        ids->Jy[i][1][k] *= 2;
        ids->Jy[i][grd->nyn - 2][k] *= 2;
        // Jz
        ids->Jz[i][1][k] *= 2;
        ids->Jz[i][grd->nyn - 2][k] *= 2;

        // pxx
        ids->pxx[i][1][k] *= 2;
        ids->pxx[i][grd->nyn - 2][k] *= 2;
        // pxy
        ids->pxy[i][1][k] *= 2;
        ids->pxy[i][grd->nyn - 2][k] *= 2;
        // pxz
        ids->pxz[i][1][k] *= 2;
        ids->pxz[i][grd->nyn - 2][k] *= 2;

        // pyy
        ids->pyy[i][1][k] *= 2;
        ids->pyy[i][grd->nyn - 2][k] *= 2;
        // pyz
        ids->pyz[i][1][k] *= 2;
        ids->pyz[i][grd->nyn - 2][k] *= 2;
        // pzz
        ids->pzz[i][1][k] *= 2;
        ids->pzz[i][grd->nyn - 2][k] *= 2;

      }  // end of loop

  }  // end of non PERIODICY

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == false) {
    for (int i = 1; i < grd->nxn - 1; i++)
      for (int j = 1; j < grd->nyn - 1; j++) {
        // rhon
        ids->rhon[i][j][1] *= 2;
        ids->rhon[i][j][grd->nzn - 2] *= 2;

        // Jx
        ids->Jx[i][j][1] *= 2;
        ids->Jx[i][j][grd->nzn - 2] *= 2;
        // Jy
        ids->Jy[i][j][1] *= 2;
        ids->Jy[i][j][grd->nzn - 2] *= 2;
        // Jz
        ids->Jz[i][j][1] *= 2;
        ids->Jz[i][j][grd->nzn - 2] *= 2;

        // pxx
        ids->pxx[i][j][1] *= 2;
        ids->pxx[i][j][grd->nzn - 2] *= 2;
        // pxy
        ids->pxy[i][j][1] *= 2;
        ids->pxy[i][j][grd->nzn - 2] *= 2;
        // pxz
        ids->pxz[i][j][1] *= 2;
        ids->pxz[i][j][grd->nzn - 2] *= 2;
        // pyy
        ids->pyy[i][j][1] *= 2;
        ids->pyy[i][j][grd->nzn - 2] *= 2;
        // pyz
        ids->pyz[i][j][1] *= 2;
        ids->pyz[i][j][grd->nzn - 2] *= 2;
        // pzz
        ids->pzz[i][j][1] *= 2;
        ids->pzz[i][j][grd->nzn - 2] *= 2;
      }
  }  // end of non X periodic
}

//////////
// POPULATE GHOST CELL ON CELL CENTERS
//////////

/** Apply BC to scalar interp quantity defined on center- Interpolation quantity
 */
void applyBCscalarDensC(FPinterp ***scalarC, grid *grd, parameters *param) {
  ///////////////////////
  ///
  ///    FACE
  ///

  // X direction
  if (param->PERIODICX == true) {
    for (int j = 1; j < grd->nyc - 1; j++)
      for (int k = 1; k < grd->nzc - 1; k++) {
        scalarC[0][j][k] = scalarC[grd->nxc - 2][j][k];
        scalarC[grd->nxc - 1][j][k] = scalarC[1][j][k];
      }
  }

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == true) {
    for (int i = 1; i < grd->nxc - 1; i++)
      for (int k = 1; k < grd->nzc - 1; k++) {
        // rhon
        scalarC[i][0][k] = scalarC[i][grd->nyc - 2][k];
        scalarC[i][grd->nyc - 1][k] = scalarC[i][1][k];
      }
  }

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == true) {
    for (int i = 1; i < grd->nxc - 1; i++)
      for (int j = 1; j < grd->nyc - 1; j++) {
        // rhon
        scalarC[i][j][0] = scalarC[i][j][grd->nzc - 2];
        scalarC[i][j][grd->nzc - 1] = scalarC[i][j][1];
      }
  }

  ///////////////////////
  ///
  ///    EDGES
  ///

  // X-EDGE
  if (param->PERIODICY == true || param->PERIODICZ == true) {
    for (int i = 1; i < (grd->nxc - 1); i++) {
      scalarC[i][grd->nyc - 1][grd->nzc - 1] = scalarC[i][1][1];
      scalarC[i][0][0] = scalarC[i][grd->nyc - 2][grd->nzc - 2];
      scalarC[i][0][grd->nzc - 1] = scalarC[i][grd->nyc - 2][1];
      scalarC[i][grd->nyc - 1][0] = scalarC[i][1][grd->nzc - 2];
    }
  }

  // Y-EDGE
  if (param->PERIODICX == true || param->PERIODICZ == true) {
    for (int i = 1; i < (grd->nyn - 1); i++) {
      scalarC[grd->nxc - 1][i][grd->nzc - 1] = scalarC[1][i][1];
      scalarC[0][i][0] = scalarC[grd->nxc - 2][i][grd->nzc - 2];
      scalarC[0][i][grd->nzc - 1] = scalarC[grd->nxc - 2][i][1];
      scalarC[grd->nxc - 1][i][0] = scalarC[1][i][grd->nzc - 2];
    }
  }

  // Z-EDGE
  if (param->PERIODICX == true || param->PERIODICY == true) {
    for (int i = 1; i < (grd->nzc - 1); i++) {
      scalarC[grd->nxc - 1][grd->nyc - 1][i] = scalarC[1][1][i];
      scalarC[0][0][i] = scalarC[grd->nxc - 2][grd->nyc - 2][i];
      scalarC[grd->nxc - 1][0][i] = scalarC[1][grd->nyc - 2][i];
      scalarC[0][grd->nyc - 1][i] = scalarC[grd->nxc - 2][1][i];
    }
  }

  // Corners
  if (param->PERIODICX == true || param->PERIODICY == true ||
      param->PERIODICZ == true) {
    scalarC[grd->nxc - 1][grd->nyc - 1][grd->nzc - 1] = scalarC[1][1][1];
    scalarC[0][grd->nyc - 1][grd->nzc - 1] = scalarC[grd->nxc - 2][1][1];
    scalarC[grd->nxc - 1][0][grd->nzc - 1] = scalarC[1][grd->nyc - 2][1];
    scalarC[0][0][grd->nzc - 1] = scalarC[grd->nxc - 2][grd->nyc - 2][1];
    scalarC[grd->nxc - 1][grd->nyc - 1][0] = scalarC[1][1][grd->nzc - 2];
    scalarC[0][grd->nyc - 1][0] = scalarC[grd->nxc - 2][1][grd->nzc - 2];
    scalarC[grd->nxc - 1][0][0] = scalarC[1][grd->nyc - 2][grd->nzc - 2];
    scalarC[0][0][0] = scalarC[grd->nxc - 2][grd->nyc - 2][grd->nzc - 2];
  }

  // FACE NON PERIODIC: PUT Neuman condition absence of something else
  // X direction
  if (param->PERIODICX == false) {
    for (int j = 0; j < grd->nyc; j++)
      for (int k = 0; k < grd->nzc; k++) {
        scalarC[0][j][k] = scalarC[1][j][k];
        scalarC[grd->nxc - 1][j][k] = scalarC[grd->nxc - 2][j][k];
      }
  }

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == false) {
    for (int i = 0; i < grd->nxc; i++)
      for (int k = 0; k < grd->nzc; k++) {
        scalarC[i][0][k] = scalarC[i][1][k];
        scalarC[i][grd->nyc - 1][k] = scalarC[i][grd->nyc - 2][k];
      }
  }

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == false) {
    for (int i = 0; i < grd->nxc; i++)
      for (int j = 0; j < grd->nyc; j++) {
        scalarC[i][j][0] = scalarC[i][j][1];
        scalarC[i][j][grd->nzc - 1] = scalarC[i][j][grd->nzc - 2];
      }
  }
}

/** Apply BC to scalar field quantity defined on center - Interpolation quantity
 */
void applyBCscalarFieldC(FPfield ***scalarC, grid *grd, parameters *param) {
  ///////////////////////
  ///
  ///    FACE
  ///

  // X direction
  if (param->PERIODICX == true) {
    for (int j = 1; j < grd->nyc - 1; j++)
      for (int k = 1; k < grd->nzc - 1; k++) {
        scalarC[0][j][k] = scalarC[grd->nxc - 2][j][k];
        scalarC[grd->nxc - 1][j][k] = scalarC[1][j][k];
      }
  }

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == true) {
    for (int i = 1; i < grd->nxc - 1; i++)
      for (int k = 1; k < grd->nzc - 1; k++) {
        // rhon
        scalarC[i][0][k] = scalarC[i][grd->nyc - 2][k];
        scalarC[i][grd->nyc - 1][k] = scalarC[i][1][k];
      }
  }

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == true) {
    for (int i = 1; i < grd->nxc - 1; i++)
      for (int j = 1; j < grd->nyc - 1; j++) {
        // rhon
        scalarC[i][j][0] = scalarC[i][j][grd->nzc - 2];
        scalarC[i][j][grd->nzc - 1] = scalarC[i][j][1];
      }
  }

  ///////////////////////
  ///
  ///    EDGES
  ///

  // X-EDGE
  if (param->PERIODICY == true || param->PERIODICZ == true) {
    for (int i = 1; i < (grd->nxc - 1); i++) {
      scalarC[i][grd->nyc - 1][grd->nzc - 1] = scalarC[i][1][1];
      scalarC[i][0][0] = scalarC[i][grd->nyc - 2][grd->nzc - 2];
      scalarC[i][0][grd->nzc - 1] = scalarC[i][grd->nyc - 2][1];
      scalarC[i][grd->nyc - 1][0] = scalarC[i][1][grd->nzc - 2];
    }
  }

  // Y-EDGE
  if (param->PERIODICX == true || param->PERIODICZ == true) {
    for (int i = 1; i < (grd->nyn - 1); i++) {
      scalarC[grd->nxc - 1][i][grd->nzc - 1] = scalarC[1][i][1];
      scalarC[0][i][0] = scalarC[grd->nxc - 2][i][grd->nzc - 2];
      scalarC[0][i][grd->nzc - 1] = scalarC[grd->nxc - 2][i][1];
      scalarC[grd->nxc - 1][i][0] = scalarC[1][i][grd->nzc - 2];
    }
  }

  // Z-EDGE
  if (param->PERIODICX == true || param->PERIODICY == true) {
    for (int i = 1; i < (grd->nzc - 1); i++) {
      scalarC[grd->nxc - 1][grd->nyc - 1][i] = scalarC[1][1][i];
      scalarC[0][0][i] = scalarC[grd->nxc - 2][grd->nyc - 2][i];
      scalarC[grd->nxc - 1][0][i] = scalarC[1][grd->nyc - 2][i];
      scalarC[0][grd->nyc - 1][i] = scalarC[grd->nxc - 2][1][i];
    }
  }

  // Corners
  if (param->PERIODICX == true || param->PERIODICY == true ||
      param->PERIODICZ == true) {
    scalarC[grd->nxc - 1][grd->nyc - 1][grd->nzc - 1] = scalarC[1][1][1];
    scalarC[0][grd->nyc - 1][grd->nzc - 1] = scalarC[grd->nxc - 2][1][1];
    scalarC[grd->nxc - 1][0][grd->nzc - 1] = scalarC[1][grd->nyc - 2][1];
    scalarC[0][0][grd->nzc - 1] = scalarC[grd->nxc - 2][grd->nyc - 2][1];
    scalarC[grd->nxc - 1][grd->nyc - 1][0] = scalarC[1][1][grd->nzc - 2];
    scalarC[0][grd->nyc - 1][0] = scalarC[grd->nxc - 2][1][grd->nzc - 2];
    scalarC[grd->nxc - 1][0][0] = scalarC[1][grd->nyc - 2][grd->nzc - 2];
    scalarC[0][0][0] = scalarC[grd->nxc - 2][grd->nyc - 2][grd->nzc - 2];
  }

  // FACE NON PERIODIC: PUT Neuman condition absence of something else
  // X direction
  if (param->PERIODICX == false) {
    for (int j = 0; j < grd->nyc; j++)
      for (int k = 0; k < grd->nzc; k++) {
        scalarC[0][j][k] = scalarC[1][j][k];
        scalarC[grd->nxc - 1][j][k] = scalarC[grd->nxc - 2][j][k];
      }
  }

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == false) {
    for (int i = 0; i < grd->nxc; i++)
      for (int k = 0; k < grd->nzc; k++) {
        scalarC[i][0][k] = scalarC[i][1][k];
        scalarC[i][grd->nyc - 1][k] = scalarC[i][grd->nyc - 2][k];
      }
  }

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == false) {
    for (int i = 0; i < grd->nxc; i++)
      for (int j = 0; j < grd->nyc; j++) {
        scalarC[i][j][0] = scalarC[i][j][1];
        scalarC[i][j][grd->nzc - 1] = scalarC[i][j][grd->nzc - 2];
      }
  }
}

/** Apply BC to scalar field quantity defined on nodes - Interpolation quantity
 */
// set to zero ghost cell
void applyBCscalarFieldCzero(FPfield ***scalarC, grid *grd, parameters *param) {
  ///////////////////////
  ///
  ///    FACE
  ///

  // X direction
  if (param->PERIODICX == true) {
    for (int j = 1; j < grd->nyc - 1; j++)
      for (int k = 1; k < grd->nzc - 1; k++) {
        scalarC[0][j][k] = scalarC[grd->nxc - 2][j][k];
        scalarC[grd->nxc - 1][j][k] = scalarC[1][j][k];
      }
  }

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == true) {
    for (int i = 1; i < grd->nxc - 1; i++)
      for (int k = 1; k < grd->nzc - 1; k++) {
        scalarC[i][0][k] = scalarC[i][grd->nyc - 2][k];
        scalarC[i][grd->nyc - 1][k] = scalarC[i][1][k];
      }
  }

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == true) {
    for (int i = 1; i < grd->nxc - 1; i++)
      for (int j = 1; j < grd->nyc - 1; j++) {
        // rhon
        scalarC[i][j][0] = scalarC[i][j][grd->nzc - 2];
        scalarC[i][j][grd->nzc - 1] = scalarC[i][j][1];
      }
  }

  ///////////////////////
  ///
  ///    EDGES
  ///

  // X-EDGE
  if (param->PERIODICY == true || param->PERIODICZ == true) {
    for (int i = 1; i < (grd->nxc - 1); i++) {
      scalarC[i][grd->nyc - 1][grd->nzc - 1] = scalarC[i][1][1];
      scalarC[i][0][0] = scalarC[i][grd->nyc - 2][grd->nzc - 2];
      scalarC[i][0][grd->nzc - 1] = scalarC[i][grd->nyc - 2][1];
      scalarC[i][grd->nyc - 1][0] = scalarC[i][1][grd->nzc - 2];
    }
  }

  // Y-EDGE
  if (param->PERIODICX == true || param->PERIODICZ == true) {
    for (int i = 1; i < (grd->nyn - 1); i++) {
      scalarC[grd->nxc - 1][i][grd->nzc - 1] = scalarC[1][i][1];
      scalarC[0][i][0] = scalarC[grd->nxc - 2][i][grd->nzc - 2];
      scalarC[0][i][grd->nzc - 1] = scalarC[grd->nxc - 2][i][1];
      scalarC[grd->nxc - 1][i][0] = scalarC[1][i][grd->nzc - 2];
    }
  }

  // Z-EDGE
  if (param->PERIODICX == true || param->PERIODICY == true) {
    for (int i = 1; i < (grd->nzc - 1); i++) {
      scalarC[grd->nxc - 1][grd->nyc - 1][i] = scalarC[1][1][i];
      scalarC[0][0][i] = scalarC[grd->nxc - 2][grd->nyc - 2][i];
      scalarC[grd->nxc - 1][0][i] = scalarC[1][grd->nyc - 2][i];
      scalarC[0][grd->nyc - 1][i] = scalarC[grd->nxc - 2][1][i];
    }
  }

  // Corners
  if (param->PERIODICX == true || param->PERIODICY == true ||
      param->PERIODICZ == true) {
    scalarC[grd->nxc - 1][grd->nyc - 1][grd->nzc - 1] = scalarC[1][1][1];
    scalarC[0][grd->nyc - 1][grd->nzc - 1] = scalarC[grd->nxc - 2][1][1];
    scalarC[grd->nxc - 1][0][grd->nzc - 1] = scalarC[1][grd->nyc - 2][1];
    scalarC[0][0][grd->nzc - 1] = scalarC[grd->nxc - 2][grd->nyc - 2][1];
    scalarC[grd->nxc - 1][grd->nyc - 1][0] = scalarC[1][1][grd->nzc - 2];
    scalarC[0][grd->nyc - 1][0] = scalarC[grd->nxc - 2][1][grd->nzc - 2];
    scalarC[grd->nxc - 1][0][0] = scalarC[1][grd->nyc - 2][grd->nzc - 2];
    scalarC[0][0][0] = scalarC[grd->nxc - 2][grd->nyc - 2][grd->nzc - 2];
  }

  // FACE NON PERIODIC: PUT Neuman condition absence of something else
  // X direction
  if (param->PERIODICX == false) {
    for (int j = 0; j < grd->nyc; j++)
      for (int k = 0; k < grd->nzc; k++) {
        scalarC[0][j][k] = 0.0;
        scalarC[grd->nxc - 1][j][k] = 0.0;
      }
  }

  // Periodic boundary conditions in Y direction
  if (param->PERIODICY == false) {
    for (int i = 0; i < grd->nxc; i++)
      for (int k = 0; k < grd->nzc; k++) {
        scalarC[i][0][k] = 0.0;
        scalarC[i][grd->nyc - 1][k] = 0.0;
      }
  }

  // Periodic boundary conditions in Z direction
  if (param->PERIODICZ == false) {
    for (int i = 0; i < grd->nxc; i++)
      for (int j = 0; j < grd->nyc; j++) {
        scalarC[i][j][0] = 0.0;
        scalarC[i][j][grd->nzc - 1] = 0.0;
      }
  }
}
