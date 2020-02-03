#include <stdio.h>
#include "gpu/InterpDensNet_gpu.cuh"

/**
 *
 * Allocate gpu memory for InterpDensSpecies
 */
void interp_DNet_alloc_device(struct interpDensNet* idn_gpu, struct grid* grd) {
  int length = sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn;

  checkCudaErrors(cudaMalloc(&idn_gpu->rhon_flat, length));
  checkCudaErrors(cudaMalloc(
      &idn_gpu->rhoc_flat, sizeof(FPinterp) * grd->nxc * grd->nyc * grd->nzc));

  checkCudaErrors(cudaMalloc(&idn_gpu->Jx_flat, length));
  checkCudaErrors(cudaMalloc(&idn_gpu->Jy_flat, length));
  checkCudaErrors(cudaMalloc(&idn_gpu->Jz_flat, length));

  checkCudaErrors(cudaMalloc(&idn_gpu->pxx_flat, length));
  checkCudaErrors(cudaMalloc(&idn_gpu->pxy_flat, length));
  checkCudaErrors(cudaMalloc(&idn_gpu->pxz_flat, length));
  checkCudaErrors(cudaMalloc(&idn_gpu->pyy_flat, length));
  checkCudaErrors(cudaMalloc(&idn_gpu->pyz_flat, length));
  checkCudaErrors(cudaMalloc(&idn_gpu->pzz_flat, length));
}

/**
 * Copy InterpDensSpecies data from cpu to gpu.
 * Args cpu_src, gpu_dest, grid
 */
void interp_DNet_copy_to_device(struct interpDensNet* cpu_src,
                                struct interpDensNet* gpu_dest,
                                struct grid* grd) {
  int length = sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn;

  checkCudaErrors(cudaMemcpy(gpu_dest->rhon_flat, cpu_src->rhon_flat, length,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_dest->rhoc_flat, cpu_src->rhoc_flat,
                             sizeof(FPinterp) * grd->nxc * grd->nyc * grd->nzc,
                             cudaMemcpyHostToDevice));

  checkCudaErrors(cudaMemcpy(gpu_dest->Jx_flat, cpu_src->Jx_flat, length,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_dest->Jy_flat, cpu_src->Jy_flat, length,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_dest->Jz_flat, cpu_src->Jz_flat, length,
                             cudaMemcpyHostToDevice));

  checkCudaErrors(cudaMemcpy(gpu_dest->pxx_flat, cpu_src->pxx_flat, length,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_dest->pxy_flat, cpu_src->pxy_flat, length,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_dest->pxz_flat, cpu_src->pxz_flat, length,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_dest->pyy_flat, cpu_src->pyy_flat, length,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_dest->pyz_flat, cpu_src->pyz_flat, length,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_dest->pzz_flat, cpu_src->pzz_flat, length,
                             cudaMemcpyHostToDevice));
}

/**
 * Copy InterpDensSpecies data from gpu to cpu.
 * Args gpu_src, cpu_dest, grid
 */
void interp_DNet_copy_to_host(struct interpDensNet* gpu_src,
                              struct interpDensNet* cpu_dest,
                              struct grid* grd) {
  int length = sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn;

  checkCudaErrors(cudaMemcpy(cpu_dest->rhon_flat, gpu_src->rhon_flat, length,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(cpu_dest->rhoc_flat, gpu_src->rhoc_flat,
                             sizeof(FPinterp) * grd->nxc * grd->nyc * grd->nzc,
                             cudaMemcpyDeviceToHost));

  checkCudaErrors(cudaMemcpy(cpu_dest->Jx_flat, gpu_src->Jx_flat, length,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(cpu_dest->Jy_flat, gpu_src->Jy_flat, length,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(cpu_dest->Jz_flat, gpu_src->Jz_flat, length,
                             cudaMemcpyDeviceToHost));

  checkCudaErrors(cudaMemcpy(cpu_dest->pxx_flat, gpu_src->pxx_flat, length,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(cpu_dest->pxy_flat, gpu_src->pxy_flat, length,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(cpu_dest->pxz_flat, gpu_src->pxz_flat, length,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(cpu_dest->pyy_flat, gpu_src->pyy_flat, length,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(cpu_dest->pyz_flat, gpu_src->pyz_flat, length,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(cpu_dest->pzz_flat, gpu_src->pzz_flat, length,
                             cudaMemcpyDeviceToHost));
}

/**
 * Dealloc interpDensSpecies on gpu
 */
void interp_DNet_dealloc_device(struct interpDensNet* idn,
                                struct interpDensNet* idn_ptr) {
  checkCudaErrors(cudaFree(idn->Jx_flat));
  checkCudaErrors(cudaFree(idn->Jy_flat));
  checkCudaErrors(cudaFree(idn->Jz_flat));

  checkCudaErrors(cudaFree(idn->pxx_flat));
  checkCudaErrors(cudaFree(idn->pxy_flat));
  checkCudaErrors(cudaFree(idn->pxz_flat));
  checkCudaErrors(cudaFree(idn->pyy_flat));
  checkCudaErrors(cudaFree(idn->pyz_flat));
  checkCudaErrors(cudaFree(idn->pzz_flat));

  checkCudaErrors(cudaFree(idn_ptr));
}

/**
 * Allocate InterpDensSpecies memory  on gpu and copy data
 */
void interp_DNet_alloc_and_copy_to_device(struct interpDensNet* idn_cpu,
                                          struct interpDensNet* idn_gpu,
                                          struct interpDensNet** idn_gpu_ptr,
                                          struct grid* grd) {
  std::memcpy(idn_gpu, idn_cpu, sizeof(interpDensNet));

  // alloc memory for array and copy
  interp_DNet_alloc_device(idn_gpu, grd);
  interp_DNet_copy_to_device(idn_cpu, idn_gpu, grd);

  // alloc memory for struct and copy
  checkCudaErrors(cudaMalloc(idn_gpu_ptr, sizeof(interpDensNet)));
  checkCudaErrors(cudaMemcpy(*idn_gpu_ptr, idn_gpu, sizeof(interpDensNet),
                             cudaMemcpyHostToDevice));
}

void setZeroSpeciesDensities_gpu(cudaStream_t* stream, struct grid* grd,
                                 struct grid* grd_gpu_ptr,
                                 struct interpDensSpecies* ids_gpu_ptr,
                                 int id) {
  checkCudaErrors(cudaMemsetAsync(
      ids_gpu_ptr->Jx_flat, 0,
      sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn, *stream));
  checkCudaErrors(cudaMemsetAsync(
      ids_gpu_ptr->Jy_flat, 0,
      sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn, *stream));
  checkCudaErrors(cudaMemsetAsync(
      ids_gpu_ptr->Jz_flat, 0,
      sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn, *stream));
  checkCudaErrors(cudaMemsetAsync(
      ids_gpu_ptr->pxx_flat, 0,
      sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn, *stream));
  checkCudaErrors(cudaMemsetAsync(
      ids_gpu_ptr->pxy_flat, 0,
      sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn, *stream));
  checkCudaErrors(cudaMemsetAsync(
      ids_gpu_ptr->pxz_flat, 0,
      sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn, *stream));
  checkCudaErrors(cudaMemsetAsync(
      ids_gpu_ptr->pyy_flat, 0,
      sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn, *stream));
  checkCudaErrors(cudaMemsetAsync(
      ids_gpu_ptr->pyz_flat, 0,
      sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn, *stream));
  checkCudaErrors(cudaMemsetAsync(
      ids_gpu_ptr->pzz_flat, 0,
      sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn, *stream));
  checkCudaErrors(cudaMemsetAsync(
      ids_gpu_ptr->rhon_flat, 0,
      sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn, *stream));

  checkCudaErrors(cudaMemsetAsync(
      ids_gpu_ptr->rhoc_flat, 0,
      sizeof(FPinterp) * grd->nxc * grd->nyc * grd->nzc, *stream));
}

/** set all the densities to zero */
__global__ void set_zero_species_densities_nodes(struct interpDensSpecies* ids,
                                                 struct grid* grd, int n) {
  int local_index = blockIdx.x * blockDim.x + threadIdx.x;
  int x = local_index % grd->nxn;
  local_index = local_index / grd->nxn;

  int y = local_index % grd->nyn;
  int z = local_index / grd->nyn;

  if (grd->nxn <= x) {
    return;
  }
  if (grd->nyn <= y) {
    return;
  }
  if (grd->nzn <= z) {
    return;
  }

  ids->rhon_flat[get_idx(x, y, z, grd->nyn, grd->nzn)] = 0.0;
  // current
  ids->Jx_flat[get_idx(x, y, z, grd->nyn, grd->nzn)] = 0.0;
  ids->Jy_flat[get_idx(x, y, z, grd->nyn, grd->nzn)] = 0.0;
  ids->Jz_flat[get_idx(x, y, z, grd->nyn, grd->nzn)] = 0.0;
  // pressure
  ids->pxx_flat[get_idx(x, y, z, grd->nyn, grd->nzn)] = 0.0;
  ids->pxy_flat[get_idx(x, y, z, grd->nyn, grd->nzn)] = 0.0;
  ids->pxz_flat[get_idx(x, y, z, grd->nyn, grd->nzn)] = 0.0;
  ids->pyy_flat[get_idx(x, y, z, grd->nyn, grd->nzn)] = 0.0;
  ids->pyz_flat[get_idx(x, y, z, grd->nyn, grd->nzn)] = 0.0;
  ids->pzz_flat[get_idx(x, y, z, grd->nyn, grd->nzn)] = 0.0;

  if (grd->nxc <= x) {
    return;
  }
  if (grd->nyc <= y) {
    return;
  }
  if (grd->nzc <= z) {
    return;
  }

  ids->rhoc_flat[get_idx(x, y, z, grd->nyc, grd->nzc)] = 0.0;
}

/** set all the densities to zero */

void sumOverSpecies_gpu(struct interpDensNet* idn_gpu_ptr,
                        struct interpDensSpecies* ids_gpu_ptr,
                        struct grid* grd_gpu_ptr, struct grid* grd,
                        struct parameters* param) {
  int blocks = (grd->nxn * grd->nyn * grd->nzn + param->threads_per_block - 1) /
               param->threads_per_block;

  sum_over_species_gpu<<<blocks, param->threads_per_block>>>(idn_gpu_ptr, ids_gpu_ptr,
                                                      grd_gpu_ptr);
}

__global__ void sum_over_species_gpu(struct interpDensNet* idn,
                                     struct interpDensSpecies* ids,
                                     struct grid* grd) {
  int local_index = blockIdx.x * blockDim.x + threadIdx.x;
  int x = local_index % grd->nxn;
  local_index = local_index / grd->nxn;

  int y = local_index % grd->nyn;
  int z = local_index / grd->nyn;

  if (grd->nxn <= x) {
    return;
  }
  if (grd->nyn <= y) {
    return;
  }
  if (grd->nzn <= z) {
    return;
  }

  // density
  atomicAdd(&idn->rhon_flat[get_idx(x, y, z, grd->nyn, grd->nzn)],
            ids->rhon_flat[get_idx(x, y, z, grd->nyn, grd->nzn)]);

  atomicAdd(&idn->Jx_flat[get_idx(x, y, z, grd->nyn, grd->nzn)],
            ids->Jx_flat[get_idx(x, y, z, grd->nyn, grd->nzn)]);
  atomicAdd(&idn->Jy_flat[get_idx(x, y, z, grd->nyn, grd->nzn)],
            ids->Jy_flat[get_idx(x, y, z, grd->nyn, grd->nzn)]);
  atomicAdd(&idn->Jz_flat[get_idx(x, y, z, grd->nyn, grd->nzn)],
            ids->Jz_flat[get_idx(x, y, z, grd->nyn, grd->nzn)]);

  atomicAdd(&idn->pxx_flat[get_idx(x, y, z, grd->nyn, grd->nzn)],
            ids->pxx_flat[get_idx(x, y, z, grd->nyn, grd->nzn)]);
  atomicAdd(&idn->pxy_flat[get_idx(x, y, z, grd->nyn, grd->nzn)],
            ids->pxy_flat[get_idx(x, y, z, grd->nyn, grd->nzn)]);
  atomicAdd(&idn->pxz_flat[get_idx(x, y, z, grd->nyn, grd->nzn)],
            ids->pxz_flat[get_idx(x, y, z, grd->nyn, grd->nzn)]);
  atomicAdd(&idn->pyy_flat[get_idx(x, y, z, grd->nyn, grd->nzn)],
            ids->pyy_flat[get_idx(x, y, z, grd->nyn, grd->nzn)]);
  atomicAdd(&idn->pyz_flat[get_idx(x, y, z, grd->nyn, grd->nzn)],
            ids->pyz_flat[get_idx(x, y, z, grd->nyn, grd->nzn)]);
  atomicAdd(&idn->pzz_flat[get_idx(x, y, z, grd->nyn, grd->nzn)],
            ids->pzz_flat[get_idx(x, y, z, grd->nyn, grd->nzn)]);
}
