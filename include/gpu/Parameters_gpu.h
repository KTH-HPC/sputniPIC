#ifndef PARAMETERS_GPU_H
#define PARAMETERS_GPU_H

#include <cuda.h>
#include <cuda_runtime.h>
#include <gpu/cuda_helper.h>

inline void param_alloc_and_copy_to_device(
    struct parameters* parameters_cpu, struct parameters** parameters_gpu_ptr) {
  checkCudaErrors(cudaMalloc(parameters_gpu_ptr, sizeof(parameters)));
  checkCudaErrors(cudaMemcpy(*parameters_gpu_ptr, parameters_cpu,
                             sizeof(parameters), cudaMemcpyHostToDevice));
};

#endif