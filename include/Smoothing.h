#ifndef SMOOTHING_H
#define SMOOTHING_H

#include "BC.h"
#include "Grid.h"
#include "Parameters.h"
#include "PrecisionTypes.h"

/** Smmoth Interpolation Quantity defined on Center */
void smoothInterpScalarC(FPinterp ***vectorC, grid *grd, parameters *param);

/** Smmoth Interpolation Quantity defined on Nodes */
void smoothInterpScalarN(FPinterp ***vectorN, grid *grd, parameters *param);

#endif
