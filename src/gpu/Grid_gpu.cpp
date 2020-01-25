
#include "gpu/Grid_gpu.h"

/**
 * Allocate memory for grid data on gpu
 */
void grid_allocate_device(struct grid* grd_gpu) {
  // allocate memory for grid point arrays
  int count = sizeof(FPfield) * grd_gpu->nxn * grd_gpu->nyn * grd_gpu->nzn;

  checkCudaErrors(cudaMalloc(&grd_gpu->XN_flat, count));
  checkCudaErrors(cudaMalloc(&grd_gpu->YN_flat, count));
  checkCudaErrors(cudaMalloc(&grd_gpu->ZN_flat, count));
}

/**
 * Copy grid data from CPU to GPU
 */
void grid_copy_to_device(struct grid* cpu_src, struct grid* gpu_dest) {
  size_t count = sizeof(FPfield) * cpu_src->nxn * cpu_src->nyn * cpu_src->nzn;

  checkCudaErrors(cudaMemcpy(gpu_dest->XN_flat, cpu_src->XN_flat, count,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_dest->YN_flat, cpu_src->YN_flat, count,
                             cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_dest->ZN_flat, cpu_src->ZN_flat, count,
                             cudaMemcpyHostToDevice));
}

/**
 * Copy grid data from CPU to GPU
 */
void grid_copy_to_host(struct grid* gpu_src, struct grid* cpu_dest) {
  int len = sizeof(FPfield) * gpu_src->nxn * gpu_src->nyn * gpu_src->nzn;

  checkCudaErrors(cudaMemcpy(cpu_dest->XN_flat, gpu_src->XN_flat, len,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(cpu_dest->YN_flat, gpu_src->YN_flat, len,
                             cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaMemcpy(cpu_dest->ZN_flat, gpu_src->ZN_flat, len,
                             cudaMemcpyDeviceToHost));
}

/** Deallocate grid*/
void grid_deallocate_device(struct grid* grid_gpu, struct grid* grid_gpu_ptr) {
  checkCudaErrors(cudaFree(grid_gpu->XN_flat));
  checkCudaErrors(cudaFree(grid_gpu->YN_flat));
  checkCudaErrors(cudaFree(grid_gpu->ZN_flat));

  checkCudaErrors(cudaFree(grid_gpu_ptr));
}

void grid_alloc_and_copy_to_device(struct grid* grid_cpu, struct grid* grid_gpu,
                                   struct grid** grid_gpu_ptr) {
  // Copy fields from CPU to GPU structs.
  std::memcpy(grid_gpu, grid_cpu, sizeof(grid));

  // Allocate arrays in the struct.
  grid_allocate_device(grid_gpu);

  // Copy data from the CPU struct to GPU struct.
  grid_copy_to_device(grid_cpu, grid_gpu);

  // Allocate the struct on the GPU.
  checkCudaErrors(cudaMalloc(grid_gpu_ptr, sizeof(grid)));

  // Copy the struct to the GPU.
  checkCudaErrors(cudaMemcpy(*grid_gpu_ptr, grid_gpu, sizeof(grid),
                             cudaMemcpyHostToDevice));
}
