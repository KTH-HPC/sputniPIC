#ifndef GRID_GPU_H
#define GRID_GPU_H

#include <cstring>
#include <iostream>

#include "Alloc.h"
#include "Grid.h"
#include "Parameters.h"
#include "PrecisionTypes.h"
#include "gpu/cuda_helper.h"

/**
 * Allocate grid data to gpu. Run only once!
 */
void grid_allocate_device(struct grid*, struct grid**);

/**
 * copy grid data from host to device
 */
void grid_copy_to_device(struct grid*, struct grid*);

/**
 * copy grid data from device to host
 */
void grid_copy_to_host(struct grid*, struct grid*);

/**
 * Deallocate grid data
 */
void grid_deallocate_device(struct grid*, struct grid*);

/**
 *
 */
void grid_alloc_and_copy_to_device(struct grid* grid_cpu, struct grid* grid_gpu,
                                   struct grid** grid_gpu_ptr);

#endif
