#include "gpu/cuda_helper.h"
#include "PrecisionTypes.h"
#include <iostream>

size_t get_appropriate_batch_size(struct parameters* param) {
  size_t free;
  size_t total;
  int num_devices;
  checkCudaErrors(cudaMemGetInfo(&free, &total));
  checkCudaErrors(cudaGetDeviceCount(&num_devices));
  std::cout << "GPU has free memory: " << free << "x" << num_devices << std::endl;
  free = free * num_devices;

  size_t batch_size = 0;
  //size_t batch_size = (long double)free / 6.0 / sizeof(FPpart);
  // Select the largest species.
  for (int is = 0; is < param->ns; is++) {
    if (param->npMax[is] > batch_size) batch_size = param->npMax[is];
  }

  batch_size /= param->number_of_batches;

  // total memory needed is num_particles*num_species*6*sizeof(float)
  while (batch_size * sizeof(FPpart) * 6 * param->ns > (free * 9) / 10) {
    std::cout << "The wanted batchsize is too large for GPU memory." << std::endl;
    // Reduce batchsize by half while we are over the memory limit
    //batch_size /= 2;
    batch_size = (long double) batch_size * 0.8;
    if (batch_size < 128) {
      // No memory left, we cannot run this simulation at this time.
      throw "Batch_size too small, GPU out of memory.";
      return -1;
    }
  }

  return batch_size;
}
