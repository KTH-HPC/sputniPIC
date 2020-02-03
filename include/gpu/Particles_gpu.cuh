#ifndef PARTICLES_GPU_H
#define PARTICLES_GPU_H

#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include "device_launch_parameters.h"

#include "Alloc.h"
#include "EMfield.h"
#include "Grid.h"
#include "InterpDensSpecies.h"
#include "Parameters.h"
#include "Particles.h"
#include "PrecisionTypes.h"

#include "cudaDummy.cuh"
#include "gpu/EMfield_gpu.h"
#include "gpu/Grid_gpu.h"
#include "gpu/cuda_helper.h"

/**
 * This struct contains the general information about the particle species.
 */
struct particles_info_gpu {
  /** species ID: 0, 1, 2 , ... */
  int species_ID;

  /** maximum number of particles of this species on this domain. used for
   * memory allocation */
  long npmax;
  /** number of particles of this species on this domain */
  long nop;

  /** Electron and ions have different number of iterations: ions moves slower
   * than ions */
  int NiterMover;
  /** number of particle of subcycles in the mover */
  int n_sub_cycles;

  /** number of particles per cell */
  int npcel;
  /** number of particles per cell - X direction */
  int npcelx;
  /** number of particles per cell - Y direction */
  int npcely;
  /** number of particles per cell - Z direction */
  int npcelz;

  /** charge over mass ratio */
  FPpart qom;

  /* drift and thermal velocities for this species */
  FPpart u0, v0, w0;
  FPpart uth, vth, wth;
};

/**
 * This struct contains the positions for particles
 */
struct particles_positions_gpu {
  /** particle arrays: 1D arrays[npmax] */
  FPpart* x;
  FPpart* y;
  FPpart* z;
  FPpart* u;
  FPpart* v;
  FPpart* w;

  /** q must have precision of interpolated quantities: typically double. Not
   * used in mover */
  FPinterp* q;
};

/**
 * allocate memory for particle data on cpu. Completely analogous with
 * particle_allocate in particles.cpp, but uses pinned memory.
 */
void particle_allocate_host(struct parameters* param, struct particles* part,
                            int is);

/** deallocate particles on host, that have been allocated with
 * particle_allocate_host */
void particle_deallocate_host(struct particles*);

/**
 * Alloc memory for particle positions on gpu
 */
void particles_positions_alloc_device(
    struct particles_positions_gpu* part_pos,
    struct particles_positions_gpu** part_pos_ptr, size_t length);

/**
 * Dealloc memory on gpu for particle positions
 */
void particles_positions_dealloc_device(
    struct particles_positions_gpu* part_pos,
    struct particles_positions_gpu* part_pos_ptr);
/**
 *
 */
void particles_positions_copy_to_device(cudaStream_t* stream, struct particles*,
                                        struct particles_positions_gpu*, size_t,
                                        size_t, size_t);
void particles_positions_copy_to_host(cudaStream_t* stream,
                                      struct particles_positions_gpu*,
                                      struct particles*, size_t, size_t,
                                      size_t);

/**
 * Allocate memory on the device for the particle species information
 */
void particles_info_alloc_and_copy_to_device(struct particles_info_gpu*,
                                             struct particles_info_gpu**,
                                             struct particles*);

/**
 * Dealloc memory for particles species information
 */
void particles_info_dealloc_device(
    struct particles_info_gpu* part_info_gpu_ptr);

/**
 * Update all particles of a species on gpu in batches of size batchsize. Copies
 * from cpu to gpu, updates positions and moves back. Batchsize should be a
 * multiple of 512
 *
 */
int batch_update_particles(cudaStream_t* stream, struct particles* part_cpu,
                           struct particles_positions_gpu* part_gpu,
                           struct particles_positions_gpu* part_gpu_ptr,
                           struct particles_info_gpu* part_info_gpu_ptr,
                           struct EMfield* field_gpu_ptr,
                           struct grid* grd_gpu_ptr,
                           struct interpDensSpecies* ids_gpu_ptr,
                           struct parameters* param_cpu,
                           struct parameters* param_gpu_ptr, int batchsize);

/**
 * update particles and interpolate to grid
 */
__global__ void move_and_interpolate(struct particles_positions_gpu* part_pos,
                                     particles_info_gpu* part_info,
                                     struct EMfield* field, struct grid* grd,
                                     struct parameters* param,
                                     struct interpDensSpecies* ids_gpu,
                                     int num_particles);

/**
 * Update particle position and velocities on gpu
 */
__device__ void mover_PC_gpu(struct particles_positions_gpu*,
                             struct particles_info_gpu*, struct EMfield*,
                             struct grid*, struct parameters*);

/**
 * Interpolate particle -> grid on gpu
 */
__device__ void interpP2G_gpu(struct particles_positions_gpu* part_pos,
                              particles_info_gpu* part_info,
                              struct interpDensSpecies* ids_gpu,
                              struct grid* grd);

__device__ void bc_particles(FPpart* pos, FPpart* vel, bool periodic,
                             double len, int index);

__device__ void atomic_add_pressure(FPpart part_dat, FPpart weight[2][2][2],
                                    FPinterp* ids_arr, FPfield invVOL, int ix,
                                    int iy, int iz, int stride_y, int stride_z);

#endif
