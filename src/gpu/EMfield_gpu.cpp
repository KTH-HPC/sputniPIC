#include "gpu/EMfield_gpu.h"

/** allocate electric and magnetic field */
void field_allocate_device(struct grid* grd, struct EMfield* field) {
  size_t count = sizeof(FPfield) * grd->nxn * grd->nyn * grd->nzn;

  // E on nodes
  checkCudaErrors(cudaMalloc(&field->Ex_flat, count));
  checkCudaErrors(cudaMalloc(&field->Ey_flat, count));
  checkCudaErrors(cudaMalloc(&field->Ez_flat, count));

  // B on nodes
  checkCudaErrors(cudaMalloc(&field->Bxn_flat, count));
  checkCudaErrors(cudaMalloc(&field->Byn_flat, count));
  checkCudaErrors(cudaMalloc(&field->Bzn_flat, count));
}

/** deallocate electric and magnetic field */
void field_deallocate_device(struct EMfield* gpu_field,
                             struct EMfield* gpu_field_ptr) {
  // Start with deallocating pointers inside the struct
  // Deallocate E-fields
  checkCudaErrors(cudaFree(gpu_field->Ex_flat));
  checkCudaErrors(cudaFree(gpu_field->Ey_flat));
  checkCudaErrors(cudaFree(gpu_field->Ez_flat));

  // Deallocate B-fields
  checkCudaErrors(cudaFree(gpu_field->Bxn_flat));
  checkCudaErrors(cudaFree(gpu_field->Byn_flat));
  checkCudaErrors(cudaFree(gpu_field->Bzn_flat));

  // Free the object.
  checkCudaErrors(cudaFree(gpu_field_ptr));
}

void copy_from_host_to_device(struct EMfield* gpu_field,
                              struct EMfield* cpu_field, size_t count) {
  checkCudaErrors(cudaMemcpy(gpu_field->Ex_flat, cpu_field->Ex_flat,
                             sizeof(FPfield) * count, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_field->Ey_flat, cpu_field->Ey_flat,
                             sizeof(FPfield) * count, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_field->Ez_flat, cpu_field->Ez_flat,
                             sizeof(FPfield) * count, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_field->Bxn_flat, cpu_field->Bxn_flat,
                             sizeof(FPfield) * count, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_field->Byn_flat, cpu_field->Byn_flat,
                             sizeof(FPfield) * count, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(gpu_field->Bzn_flat, cpu_field->Bzn_flat,
                             sizeof(FPfield) * count, cudaMemcpyHostToDevice));
}

void field_alloc_and_copy_to_device(struct grid* grid,
                                    struct EMfield* gpu_field,
                                    struct EMfield* cpu_field,
                                    struct EMfield** gpu_field_ptr) {
  // Allocates the pointers inside the struct on the GPU.
  field_allocate_device(grid, gpu_field);
  size_t count = grid->nxn * grid->nyn * grid->nzn;

  // Copy data to the pointers in the struct.
  copy_from_host_to_device(gpu_field, cpu_field, count);

  // Allocate the struct on the GPU.
  checkCudaErrors(cudaMalloc(gpu_field_ptr, sizeof(EMfield)));

  // Copy the struct to the GPU.
  checkCudaErrors(cudaMemcpy(*gpu_field_ptr, gpu_field, sizeof(EMfield),
                             cudaMemcpyHostToDevice));
}
