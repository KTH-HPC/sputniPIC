#include "gpu/cuda_helper.h"
#include "PrecisionTypes.h"


int get_appropriate_batch_size(int ns)
{
    size_t free;
    size_t total;
    checkCudaErrors(cudaMemGetInfo (&free, &total));

    //int batch_size = (long double)free / 6.0 / sizeof(FPpart);
    int batch_size = THREADS_PER_BLOCK * 4096;
    // total memory needed is num_particles*num_species*6*sizeof(float)
    while(batch_size * sizeof(FPpart) * 6 * ns > (free * 9)/10){

        // Reduce batchsize by half while we are over the memory limit
        batch_size /= 2;
        if(batch_size < 128){
            // No memory left, we cannot run this simulation at this time.
            throw "Batch_size too small, GPU out of memory.";
            return -1;
        }
    }
    return batch_size;
}
