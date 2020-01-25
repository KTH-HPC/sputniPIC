#ifndef INTERP_DENS_SPECIES_GPU_H
#define INTERP_DENS_SPECIES_GPU_H

#include "InterpDensSpecies.h"

/**
 *
 * Allocate gpu memory for InterpDensSpecies
 */
void interp_DS_alloc_device(struct interpDensSpecies* ids_gpu,
                            struct grid* grd);

/**
 * Copy InterpDensSpecies data from cpu to gpu.
 * Args cpu_src, gpu_dest, grid
 */
void interp_DS_copy_to_device(struct interpDensSpecies* cpu_src,
                              struct interpDensSpecies* gpu_dest,
                              struct grid* grd);

/**
 * Copy InterpDensSpecies data from gpu to cpu.
 * Args gpu_src, cpu_dest, grid
 */
void interp_DS_copy_to_host(struct interpDensSpecies* gpu_src,
                            struct interpDensSpecies* cpu_dest,
                            struct grid* grd);

/**
 * Dealloc interpDensSpecies on gpu
 */
void interp_DS_dealloc_device(struct interpDensSpecies* ids,
                              struct interpDensSpecies* ids_ptr);

/**
 * Allocate InterpDensSpecies memory  on gpu and copy data
 */
void interp_DS_alloc_and_copy_to_device(struct interpDensSpecies* ids_cpu,
                                        struct interpDensSpecies* ids_gpu,
                                        struct interpDensSpecies** ids_gpu_ptr,
                                        struct grid* grd);

#endif