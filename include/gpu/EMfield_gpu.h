#ifndef EMFIELD_GPU_H
#define EMFIELD_GPU_H

#include "Alloc.h"
#include "EMfield.h"
#include "Grid.h"
#include "gpu/cuda_helper.h"

/** allocate electric and magnetic field on gpu */
void field_allocate_device(struct grid*, struct EMfield*);

/** deallocate electric and magnetic field on gpu */
void field_deallocate_device(struct EMfield* gpu_field,
                             struct EMfield* gpu_field_ptr);

/** Copy field data from CPU variant to GPU variant */
void copy_from_host_to_device(struct EMfield* gpu_field,
                              struct EMfield* cpu_field, size_t count);

/** Field allocation and copy to GPU */
void field_alloc_and_copy_to_device(struct grid* grid,
                                    struct EMfield* gpu_field,
                                    struct EMfield* cpu_field,
                                    struct EMfield** gpu_field_ptr);

#endif