#ifndef RW_IO_H
#define RW_IO_H

#include "EMfield.h"
#include "Grid.h"
#include "InterpDensNet.h"
#include "InterpDensSpecies.h"
#include "Parameters.h"

/** read the inputfile given via the command line */
void readInputFile(struct parameters *param, 
                    struct directories *paths,
                    int argc, char **argv);

/** Print Simulation Parameters */
void printParameters(struct parameters *param, struct directories *paths);

/** Save Simulation Parameters */
void saveParameters(struct parameters *param, struct directories *paths);

void VTK_Write_Vectors(int cycle, struct grid *grd, struct EMfield *field, struct directories *paths);

void VTK_Write_Scalars(int cycle, struct grid *grd,
                       struct interpDensSpecies *ids,
                       struct interpDensNet *idn,
                       struct directories *paths);

#endif
