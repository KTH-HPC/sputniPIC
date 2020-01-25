#ifndef BC_GPU_H
#define BC_GPU_H

#include "Grid.h"
#include "InterpDensSpecies.h"
#include "Parameters.h"
#include "PrecisionTypes.h"

#include <cuda_runtime.h>
#include "cudaDummy.cuh"
#include "gpu/InterpDensSpecies_gpu.h"
#include "gpu/cuda_helper.h"

/**
 * Apply boundary conditions on gpu
 */
void applyBCids_gpu(cudaStream_t* stream, struct interpDensSpecies* ids_gpu,
                    struct grid* grd_gpu, struct grid* grd,
                    struct parameters* param);

__global__ void applyBCids_x_p_gpu(struct interpDensSpecies* ids,
                                   struct grid* grd);

__global__ void applyBCids_y_p_gpu(struct interpDensSpecies* ids,
                                   struct grid* grd);

__global__ void applyBCids_z_p_gpu(struct interpDensSpecies* ids,
                                   struct grid* grd);

__global__ void applyBCids_x_np_gpu(struct interpDensSpecies* ids,
                                    struct grid* grd);

__global__ void applyBCids_y_np_gpu(struct interpDensSpecies* ids,
                                    struct grid* grd);

__global__ void applyBCids_z_np_gpu(struct interpDensSpecies* ids,
                                    struct grid* grd);

__host__ __device__ void apply_periodic_bc(struct interpDensSpecies* ids,
                                           int x1, int y1, int z1, int x2,
                                           int y2, int z2, int stride_y,
                                           int stride_z);

__host__ __device__ void apply_nonperiodic_bc(struct interpDensSpecies* ids,
                                              int x1, int y1, int z1, int x2,
                                              int y2, int z2, int stride_y,
                                              int stride_z);

__host__ __device__ inline void atomic_add_and_update(FPpart* array, int x1,
                                                      int y1, int z1, int x2,
                                                      int y2, int z2,
                                                      int stride_y,
                                                      int stride_z);

__host__ __device__ inline void atomic_double_update(FPpart* array, int x1,
                                                     int y1, int z1, int x2,
                                                     int y2, int z2,
                                                     int stride_y,
                                                     int stride_z);

#endif