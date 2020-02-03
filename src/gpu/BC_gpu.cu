#include "gpu/BC_gpu.cuh"
//#include<stdio.h>

void applyBCids_gpu(cudaStream_t* stream, struct interpDensSpecies* ids_gpu,
                    struct grid* grd_gpu, struct grid* grd,
                    struct parameters* param) {
  int th_x = grd->nyn * grd->nzn;
  int th_y = grd->nxn * grd->nzn;
  int th_z = grd->nxn * grd->nyn;

  // printf("thx:%d, thy:%d, thz:%d\n", th_x, th_y, th_z);

  if (param->PERIODICX) {
    applyBCids_x_p_gpu<<<(th_x + param->threads_per_block - 1) / param->threads_per_block,
                         param->threads_per_block>>>(ids_gpu, grd_gpu);
  } else {
    applyBCids_x_np_gpu<<<(th_x + param->threads_per_block - 1) / param->threads_per_block,
                          param->threads_per_block>>>(ids_gpu, grd_gpu);
  }
  if (param->PERIODICY) {
    applyBCids_y_p_gpu<<<(th_y + param->threads_per_block - 1) / param->threads_per_block,
                         param->threads_per_block>>>(ids_gpu, grd_gpu);
  } else {
    applyBCids_y_np_gpu<<<(th_y + param->threads_per_block - 1) / param->threads_per_block,
                          param->threads_per_block>>>(ids_gpu, grd_gpu);
  }
  if (param->PERIODICZ) {
    applyBCids_z_p_gpu<<<(th_z + param->threads_per_block - 1) / param->threads_per_block,
                         param->threads_per_block>>>(ids_gpu, grd_gpu);
  } else {
    applyBCids_z_np_gpu<<<(th_z + param->threads_per_block - 1) / param->threads_per_block,
                          param->threads_per_block>>>(ids_gpu, grd_gpu);
  }
}

__global__ void applyBCids_x_p_gpu(struct interpDensSpecies* ids,
                                   struct grid* grd) {
  int local_index = blockIdx.x * blockDim.x + threadIdx.x;

  int y = local_index % grd->nyn;
  int z = local_index / grd->nyn;

  if (y < 1 || grd->nyn - 1 <= y) {
    return;
  }
  if (z < 1 || grd->nzn - 1 <= z) {
    return;
  }

  apply_periodic_bc(ids, 1, y, z, grd->nxn - 2, y, z, grd->nyn, grd->nzn);
}

__global__ void applyBCids_y_p_gpu(struct interpDensSpecies* ids,
                                   struct grid* grd) {
  int local_index = blockIdx.x * blockDim.x + threadIdx.x;

  int x = local_index % grd->nxn;
  int z = local_index / grd->nxn;

  if (x < 1 || grd->nxn - 1 <= x) {
    return;
  }
  if (z < 1 || grd->nzn - 1 <= z) {
    return;
  }

  apply_periodic_bc(ids, x, 1, z, x, grd->nyn - 2, z, grd->nyn, grd->nzn);
}

__global__ void applyBCids_z_p_gpu(struct interpDensSpecies* ids,
                                   struct grid* grd) {
  int local_index = blockIdx.x * blockDim.x + threadIdx.x;

  int x = local_index % grd->nxn;
  int y = local_index / grd->nxn;

  if (x < 1 || grd->nxn - 1 <= x) {
    return;
  }
  if (y < 1 || grd->nyn - 1 <= y) {
    return;
  }

  apply_periodic_bc(ids, x, y, 1, x, y, grd->nzn - 2, grd->nyn, grd->nzn);
}

__global__ void applyBCids_x_np_gpu(struct interpDensSpecies* ids,
                                    struct grid* grd) {
  int local_index = blockIdx.x * blockDim.x + threadIdx.x;

  int y = local_index % grd->nyn;
  int z = local_index / grd->nyn;

  if (y < 1 || grd->nyn - 1 <= y) {
    return;
  }
  if (z < 1 || grd->nzn - 1 <= z) {
    return;
  }

  apply_nonperiodic_bc(ids, 1, y, z, grd->nxn - 2, y, z, grd->nyn, grd->nzn);
}

__global__ void applyBCids_y_np_gpu(struct interpDensSpecies* ids,
                                    struct grid* grd) {
  int local_index = blockIdx.x * blockDim.x + threadIdx.x;

  int x = local_index % grd->nxn;
  int z = local_index / grd->nxn;

  if (x < 1 || grd->nxn - 1 <= x) {
    return;
  }
  if (z < 1 || grd->nzn - 1 <= z) {
    return;
  }

  apply_nonperiodic_bc(ids, x, 1, z, x, grd->nyn - 2, z, grd->nyn, grd->nzn);
}

__global__ void applyBCids_z_np_gpu(struct interpDensSpecies* ids,
                                    struct grid* grd) {
  int local_index = blockIdx.x * blockDim.x + threadIdx.x;

  int x = local_index % grd->nxn;
  int y = local_index / grd->nxn;

  if (x < 1 || grd->nxn - 1 <= x) {
    return;
  }
  if (y < 1 || grd->nyn - 1 <= y) {
    return;
  }
  apply_nonperiodic_bc(ids, x, y, 1, x, y, grd->nzn - 2, grd->nyn, grd->nzn);
}

__host__ __device__ void apply_periodic_bc(struct interpDensSpecies* ids,
                                           int x1, int y1, int z1, int x2,
                                           int y2, int z2, int stride_y,
                                           int stride_z) {
  atomic_add_and_update(ids->rhon_flat, x1, y1, z1, x2, y2, z2, stride_y,
                        stride_z);

  atomic_add_and_update(ids->Jx_flat, x1, y1, z1, x2, y2, z2, stride_y,
                        stride_z);
  atomic_add_and_update(ids->Jy_flat, x1, y1, z1, x2, y2, z2, stride_y,
                        stride_z);
  atomic_add_and_update(ids->Jz_flat, x1, y1, z1, x2, y2, z2, stride_y,
                        stride_z);

  atomic_add_and_update(ids->pxx_flat, x1, y1, z1, x2, y2, z2, stride_y,
                        stride_z);
  atomic_add_and_update(ids->pyy_flat, x1, y1, z1, x2, y2, z2, stride_y,
                        stride_z);
  atomic_add_and_update(ids->pzz_flat, x1, y1, z1, x2, y2, z2, stride_y,
                        stride_z);

  atomic_add_and_update(ids->pxy_flat, x1, y1, z1, x2, y2, z2, stride_y,
                        stride_z);
  atomic_add_and_update(ids->pyz_flat, x1, y1, z1, x2, y2, z2, stride_y,
                        stride_z);
  atomic_add_and_update(ids->pxz_flat, x1, y1, z1, x2, y2, z2, stride_y,
                        stride_z);
}

__host__ __device__ void apply_nonperiodic_bc(struct interpDensSpecies* ids,
                                              int x1, int y1, int z1, int x2,
                                              int y2, int z2, int stride_y,
                                              int stride_z) {
  atomic_double_update(ids->rhon_flat, x1, y1, z1, x2, y2, z2, stride_y,
                       stride_z);

  atomic_double_update(ids->Jx_flat, x1, y1, z1, x2, y2, z2, stride_y,
                       stride_z);
  atomic_double_update(ids->Jy_flat, x1, y1, z1, x2, y2, z2, stride_y,
                       stride_z);
  atomic_double_update(ids->Jz_flat, x1, y1, z1, x2, y2, z2, stride_y,
                       stride_z);

  atomic_double_update(ids->pxx_flat, x1, y1, z1, x2, y2, z2, stride_y,
                       stride_z);
  atomic_double_update(ids->pyy_flat, x1, y1, z1, x2, y2, z2, stride_y,
                       stride_z);
  atomic_double_update(ids->pzz_flat, x1, y1, z1, x2, y2, z2, stride_y,
                       stride_z);

  atomic_double_update(ids->pxy_flat, x1, y1, z1, x2, y2, z2, stride_y,
                       stride_z);
  atomic_double_update(ids->pyz_flat, x1, y1, z1, x2, y2, z2, stride_y,
                       stride_z);
  atomic_double_update(ids->pxz_flat, x1, y1, z1, x2, y2, z2, stride_y,
                       stride_z);
}

__host__ __device__ inline void atomic_add_and_update(FPpart* array, int x1,
                                                      int y1, int z1, int x2,
                                                      int y2, int z2,
                                                      int stride_y,
                                                      int stride_z) {
  // atomicAdd(&array[get_idx(x1, y1, z1, stride_y, stride_z)],
  // array[get_idx(x2, y2, z2, stride_y, stride_z)]);
  // atomicExch(&array[get_idx(x2, y2, z2, stride_y, stride_z)],
  // array[get_idx(x1, y1, z1, stride_y, stride_z)]);

  array[get_idx(x1, y1, z1, stride_y, stride_z)] +=
      array[get_idx(x2, y2, z2, stride_y, stride_z)];
  array[get_idx(x2, y2, z2, stride_y, stride_z)] =
      array[get_idx(x1, y1, z1, stride_y, stride_z)];
}

__host__ __device__ inline void atomic_double_update(FPpart* array, int x1,
                                                     int y1, int z1, int x2,
                                                     int y2, int z2,
                                                     int stride_y,
                                                     int stride_z) {
  array[get_idx(x1, y1, z1, stride_y, stride_z)] *= 2;
  array[get_idx(x2, y2, z2, stride_y, stride_z)] *= 2;

  // atomicAdd(&array[get_idx(x1,y1,z1, stride_y, stride_z)],
  // array[get_idx(x1,y1,z1, stride_y, stride_z)]);
  // atomicAdd(&array[get_idx(x2,y2,z2, stride_y, stride_z)],
  // array[get_idx(x2,y2,z2, stride_y, stride_z)]);
}
