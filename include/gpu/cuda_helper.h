#ifndef CUDA_HELPER_H
#define CUDA_HELPER_H

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "Parameters.h"

#ifdef __INTELLISENSE__
#include "gpu/intellisence_cuda_intrinsic.h"
#endif

/* from cuda samples */
inline void checkGpuError(cudaError_t result, char const *const func,
                          const char *const file, int const line) {
  if (result != cudaSuccess) {
    fprintf(stderr, "Cuda failure %s:%d: '%s'\n", file, line,
            cudaGetErrorString(result));
    exit(1);
  }
}

// This will output the proper CUDA error strings in the event
// that a CUDA host call returns an error'
#ifdef MEMCHECK
#define checkCudaErrors(val) checkGpuError((val), #val, __FILE__, __LINE__)
#else
#define checkCudaErrors(val) (val)
#endif

/**
 * Return an appropriate size for particle batch
 */
size_t get_appropriate_batch_size(struct parameters* param);

#endif
