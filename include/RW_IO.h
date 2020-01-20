#ifndef RW_IO_H
#define RW_IO_H

#include "Parameters.h"
#include "InterpDensNet.h"
#include "InterpDensSpecies.h"
#include "Grid.h"
#include "EMfield.h"

/** read the inputfile given via the command line */
void readInputFile(struct parameters* param, int argc, char **argv);

/** Print Simulation Parameters */
void printParameters(struct parameters* param);

/** Save Simulation Parameters */
void saveParameters(struct parameters* param);

void VTK_Write_Vectors(int cycle, struct grid *grd, struct EMfield* field);

void VTK_Write_Scalars(int cycle, struct grid *grd, struct interpDensSpecies* ids, struct interpDensNet* idn);

#endif
