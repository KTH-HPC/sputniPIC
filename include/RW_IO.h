#ifndef RW_IO_H
#define RW_IO_H

#include "EMfield.h"
#include "Grid.h"
#include "InterpDensNet.h"
#include "InterpDensSpecies.h"
#include "Parameters.h"
#include "Particles.h"

#include <string>
#include <sstream>
#include <hdf5.h>

/** read the inputfile given via the command line */
void readInputFile(struct parameters *param, int argc, char **argv);

/** Print Simulation Parameters */
void printParameters(struct parameters *param);

/** Save particle positions and energy **/
void saveParticlePositions(struct parameters *param, struct particles *part, int cycle, int species); 

/** Save Simulation Parameters */
void saveParameters(struct parameters *param);

void VTK_Write_Vectors(int cycle, struct grid *grd, struct EMfield *field, struct parameters *param);

void VTK_Write_Scalars(int cycle, struct grid *grd,
                       struct interpDensSpecies *ids,
                       struct interpDensNet *idn,
                       struct parameters *param);

void VTK_Write_Vectors_Binary(int cycle, struct grid *grd, struct EMfield *field, struct parameters *param);

void VTK_Write_Scalars_Binary(int cycle, struct grid *grd,
                              struct interpDensSpecies *ids,
                              struct interpDensNet *idn,
                              struct parameters *param);

void HDF5_Write_Particles(int cycle, struct particles *part_local, struct parameters *param);

#endif
