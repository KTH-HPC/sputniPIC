#ifndef MPI_COMM_H
#define MPI_COMM_H

#include <mpi.h>
#include <stdio.h>

#include "EMfield.h"
#include "Grid.h"
#include "InterpDensNet.h"
#include "Parameters.h"
#include "Particles.h"
#include "PrecisionTypes.h"

void mpi_reduce_densities(struct grid *, struct interpDensNet *);
void mpi_reduce_densities(struct grid *, struct interpDensSpecies *);
void mpi_broadcast_field(struct grid *, struct EMfield *);
void mpi_scatter_particles(particles *part_global, particles *part_local);

#endif