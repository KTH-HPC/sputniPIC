#ifndef INTERP_DENS_NET_GPU_H
#define INTERP_DENS_NET_GPU_H

#include "Alloc.h"
#include "Grid.h"
#include "InterpDensNet.h"
#include "PrecisionTypes.h"

#include <cstring>
#include "gpu/cuda_helper.h"

/**
 *
 * Allocate gpu memory for InterpDensSpecies
 */
void interp_DNet_alloc_device(struct interpDensNet* idn_gpu, struct grid* grd);

/**
 * Copy interpDensNet data from cpu to gpu.
 * Args cpu_src, gpu_dest, grid
 */
void interp_DNet_copy_to_device(struct interpDensNet* cpu_src,
                                struct interpDensNet* gpu_dest,
                                struct grid* grd);

/**
 * Copy interpDensNet data from gpu to cpu.
 * Args gpu_src, cpu_dest, grid
 */
void interp_DNet_copy_to_host(struct interpDensNet* gpu_src,
                              struct interpDensNet* cpu_dest, struct grid* grd);

/**
 * Dealloc interpDensNet on gpu
 */
void interp_DNet_dealloc_device(struct interpDensNet* idn,
                                struct interpDensNet* idn_ptr);

/**
 * Allocate interpDensNet memory  on gpu and copy data
 */
void interp_DNet_alloc_and_copy_to_device(struct interpDensNet* idn_cpu,
                                          struct interpDensNet* idn_gpu,
                                          struct interpDensNet** idn_gpu_ptr,
                                          struct grid* grd);

void setZeroSpeciesDensities_gpu(cudaStream_t* stream, struct grid* grd,
                                 struct grid* grd_gpu_ptr,
                                 struct interpDensSpecies* ids_gpu_ptr, int id);

__global__ void set_zero_species_densities_nodes(struct interpDensSpecies* ids,
                                                 struct grid* grd, int n);

void sumOverSpecies_gpu(struct interpDensNet* idn_gpu_ptr,
                        struct interpDensSpecies* ids_gpu_ptr,
                        struct grid* grd_gpu_ptr, struct grid* grd, struct parameters* param);

__global__ void sum_over_species_gpu(struct interpDensNet* idn,
                                     struct interpDensSpecies* ids,
                                     struct grid* grd);

#endif
