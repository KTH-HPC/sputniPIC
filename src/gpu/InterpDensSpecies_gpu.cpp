#include "gpu/InterpDensSpecies_gpu.h"
#include "PrecisionTypes.h"

#include <cstring>
#include "gpu/cuda_helper.h"

/**
 *
 * Allocate gpu memory for InterpDensSpecies
 */
void interp_DS_alloc_device(struct interpDensSpecies* ids_gpu,
                            struct grid* grd) {
  int length = sizeof(FPinterp) * grd->nxn * grd->nyn * grd->nzn;

  checkCudaErrors(cudaMalloc(&ids_gpu->rhon_flat, length));
  checkCudaErrors(cudaMalloc(
      &ids_gpu->rhoc_flat, sizeof(FPinterp) * grd->nxc * grd->nyc * grd->nzc));

  checkCudaErrors(cudaMalloc(&ids_gpu->Jx_flat, length));
  checkCudaErrors(cudaMalloc(&ids_gpu->Jy_flat, length));
  checkCudaErrors(cudaMalloc(&ids_gpu->Jz_flat, length));

  checkCudaErrors(cudaMalloc(&ids_gpu->pxx_flat, length));
  checkCudaErrors(cudaMalloc(&ids_gpu->pxy_flat, length));
  checkCudaErrors(cudaMalloc(&ids_gpu->pxz_flat, length));
  checkCudaErrors(cudaMalloc(&ids_gpu->pyy_flat, length));
  checkCudaErrors(cudaMalloc(&ids_gpu->pyz_flat, length));
  checkCudaErrors(cudaMalloc(&ids_gpu->pzz_flat, length));
}

/**
 * Copy InterpDensSpecies data from cpu to gpu.
 * Args cpu_src, gpu_dest, grid
 */
void interp_DS_copy_to_device(struct interpDensSpecies* cpu_src,
                              struct interpDensSpecies* gpu_dest,
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
void interp_DS_copy_to_host(struct interpDensSpecies* gpu_src,
                            struct interpDensSpecies* cpu_dest,
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
void interp_DS_dealloc_device(struct interpDensSpecies* ids,
                              struct interpDensSpecies* ids_ptr) {
  checkCudaErrors(cudaFree(ids->Jx_flat));
  checkCudaErrors(cudaFree(ids->Jy_flat));
  checkCudaErrors(cudaFree(ids->Jz_flat));

  checkCudaErrors(cudaFree(ids->pxx_flat));
  checkCudaErrors(cudaFree(ids->pxy_flat));
  checkCudaErrors(cudaFree(ids->pxz_flat));
  checkCudaErrors(cudaFree(ids->pyy_flat));
  checkCudaErrors(cudaFree(ids->pyz_flat));
  checkCudaErrors(cudaFree(ids->pzz_flat));

  checkCudaErrors(cudaFree(ids_ptr));
}

/**
 * Allocate InterpDensSpecies memory  on gpu and copy data
 */
void interp_DS_alloc_and_copy_to_device(struct interpDensSpecies* ids_cpu,
                                        struct interpDensSpecies* ids_gpu,
                                        struct interpDensSpecies** ids_gpu_ptr,
                                        struct grid* grd) {
  std::memcpy(ids_gpu, ids_cpu, sizeof(interpDensSpecies));

  // alloc memory for array and copy
  interp_DS_alloc_device(ids_gpu, grd);
  interp_DS_copy_to_device(ids_cpu, ids_gpu, grd);

  // alloc memory for struct and copy
  checkCudaErrors(cudaMalloc(ids_gpu_ptr, sizeof(interpDensSpecies)));
  checkCudaErrors(cudaMemcpy(*ids_gpu_ptr, ids_gpu, sizeof(interpDensSpecies),
                             cudaMemcpyHostToDevice));
}
