#include "Smoothing.h"
#include "BC.h"

/** Smmoth Interpolation Quantity defined on Center */
void smoothInterpScalarC(FPinterp ***vectorC, grid *grd, parameters *param) {
  // center cells
  int nxc = grd->nxc;
  int nyc = grd->nyc;
  int nzc = grd->nzc;

  FPinterp valueS = 0.0;
  FPinterp alpha = 0.0;

  // if smoothing is on
  if (param->SmoothON) {
    // repeat Smoothing SmoothTimes
    for (int icount = 1; icount < param->SmoothTimes + 1; icount++) {
      applyBCscalarDensC(vectorC, grd,
                         param);  // set BC on ghost cells before smooth
      // alternate Smooth Values
      if (icount % 2 == 1)
        valueS = 0.0;
      else
        valueS = param->SmoothValue;

      // weights of stencil points
      alpha = (1.0 - valueS) / 6;  // Stencil of 6

      //
      for (int i = 1; i < nxc - 1; i++)
        for (int j = 1; j < nyc - 1; j++)
          for (int k = 1; k < nzc - 1; k++)
            vectorC[i][j][k] =
                valueS * vectorC[i][j][k] +
                alpha * (vectorC[i - 1][j][k] + vectorC[i + 1][j][k] +
                         vectorC[i][j - 1][k] + vectorC[i][j + 1][k] +
                         vectorC[i][j][k - 1] + vectorC[i][j][k + 1]);
    }

  }  // IF smooth is ON

  // set BC on ghost cells before leaving the function
  applyBCscalarDensC(vectorC, grd, param);
}

/** Smmoth Interpolation Quantity defined on Nodes */
void smoothInterpScalarN(FPinterp ***vectorN, grid *grd, parameters *param) {
  // nodes
  int nxn = grd->nxn;
  int nyn = grd->nyn;
  int nzn = grd->nzn;

  FPinterp valueS = 0.0;
  FPinterp alpha = 0.0;

  // if smoothing is on
  if (param->SmoothON) {
    // repeat Smoothing SmoothTimes
    for (int icount = 1; icount < param->SmoothTimes + 1; icount++) {
      applyBCscalarDensN(vectorN, grd,
                         param);  // set BC on ghost cells before smooth

      // alternate Smooth Values
      if (icount % 2 == 1)
        valueS = 0.0;
      else
        valueS = param->SmoothValue;

      // weights of stencil points
      alpha = (1.0 - valueS) / 6;  // Stencil of 6

      //
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            vectorN[i][j][k] =
                valueS * vectorN[i][j][k] +
                alpha * (vectorN[i - 1][j][k] + vectorN[i + 1][j][k] +
                         vectorN[i][j - 1][k] + vectorN[i][j + 1][k] +
                         vectorN[i][j][k - 1] + vectorN[i][j][k + 1]);
    }

  }  // IF smooth is ON

  // set BC on ghost cells before leaving the function
  applyBCscalarDensN(vectorN, grd, param);
}
