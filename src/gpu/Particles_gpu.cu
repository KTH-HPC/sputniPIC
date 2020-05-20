#include <stdio.h>
#include "PrecisionTypes.h"
#include "gpu/Particles_gpu.cuh"

/** allocate particle arrays */
void particle_allocate_host(struct parameters* param, struct particles* part,
														int is) {
	// set species ID
	part->species_ID = is;
	// number of particles
	part->nop = param->np[is];
	// maximum number of particles
	part->npmax = param->npMax[is];

	// choose a different number of mover iterations for ions and electrons
	part->NiterMover = param->NiterMover;
	part->n_sub_cycles = param->n_sub_cycles;

	// particles per cell
	part->npcelx = param->npcelx[is];
	part->npcely = param->npcely[is];
	part->npcelz = param->npcelz[is];
	part->npcel = part->npcelx * part->npcely * part->npcelz;

	// cast it to required precision
	part->qom = (FPpart)param->qom[is];

	long npmax = part->npmax;

	// initialize drift and thermal velocities
	// drift
	part->u0 = (FPpart)param->u0[is];
	part->v0 = (FPpart)param->v0[is];
	part->w0 = (FPpart)param->w0[is];
	// thermal
	part->uth = (FPpart)param->uth[is];
	part->vth = (FPpart)param->vth[is];
	part->wth = (FPpart)param->wth[is];

	//////////////////////////////
	/// ALLOCATION PARTICLE ARRAYS
	//////////////////////////////
	checkCudaErrors(cudaMallocHost(&part->x, sizeof(FPpart) * npmax));
	checkCudaErrors(cudaMallocHost(&part->y, sizeof(FPpart) * npmax));
	checkCudaErrors(cudaMallocHost(&part->z, sizeof(FPpart) * npmax));

	checkCudaErrors(cudaMallocHost(&part->u, sizeof(FPpart) * npmax));
	checkCudaErrors(cudaMallocHost(&part->v, sizeof(FPpart) * npmax));
	checkCudaErrors(cudaMallocHost(&part->w, sizeof(FPpart) * npmax));

	checkCudaErrors(cudaMallocHost(&part->q, sizeof(FPinterp) * npmax));
}

/** deallocate particles allocated with */
void particle_deallocate_host(struct particles* part) {
	checkCudaErrors(cudaFreeHost(part->x));
	checkCudaErrors(cudaFreeHost(part->y));
	checkCudaErrors(cudaFreeHost(part->z));

	checkCudaErrors(cudaFreeHost(part->u));
	checkCudaErrors(cudaFreeHost(part->v));
	checkCudaErrors(cudaFreeHost(part->w));

	checkCudaErrors(cudaFreeHost(part->q));
}

void particles_positions_alloc_device(
		struct particles_positions_gpu* part_pos,
		struct particles_positions_gpu** part_pos_ptr, size_t length) {
	checkCudaErrors(cudaMalloc(part_pos_ptr, sizeof(particles_positions_gpu)));

	checkCudaErrors(cudaMalloc(&part_pos->x, length * sizeof(FPpart)));
	checkCudaErrors(cudaMalloc(&part_pos->y, length * sizeof(FPpart)));
	checkCudaErrors(cudaMalloc(&part_pos->z, length * sizeof(FPpart)));
	checkCudaErrors(cudaMalloc(&part_pos->u, length * sizeof(FPpart)));
	checkCudaErrors(cudaMalloc(&part_pos->v, length * sizeof(FPpart)));
	checkCudaErrors(cudaMalloc(&part_pos->w, length * sizeof(FPpart)));
	checkCudaErrors(cudaMalloc(&part_pos->q, length * sizeof(FPinterp)));

	checkCudaErrors(cudaMemcpy(*part_pos_ptr, part_pos,
														 sizeof(particles_positions_gpu),
														 cudaMemcpyHostToDevice));
}

void particles_positions_dealloc_device(
		struct particles_positions_gpu* part_pos,
		struct particles_positions_gpu* part_pos_ptr) {
	checkCudaErrors(cudaFree(part_pos_ptr));

	checkCudaErrors(cudaFree(part_pos->x));
	checkCudaErrors(cudaFree(part_pos->y));
	checkCudaErrors(cudaFree(part_pos->z));
	checkCudaErrors(cudaFree(part_pos->u));
	checkCudaErrors(cudaFree(part_pos->v));
	checkCudaErrors(cudaFree(part_pos->w));
	checkCudaErrors(cudaFree(part_pos->q));
}

/**
 * Copy num_particles number of particles from cpu memory to gpu
 */
void particles_positions_copy_to_device(
		cudaStream_t* stream, 
		struct particles* particles,
		struct particles_positions_gpu* part_pos, 
		size_t idx_start_cpu,
		size_t idx_start_gpu,
		size_t num_particles
		) {

	checkCudaErrors(cudaMemcpyAsync(
			part_pos->x + idx_start_gpu, particles->x + idx_start_cpu,
			num_particles * sizeof(FPpart), cudaMemcpyHostToDevice, *stream));
	checkCudaErrors(cudaMemcpyAsync(
			part_pos->y + idx_start_gpu, particles->y + idx_start_cpu,
			num_particles * sizeof(FPpart), cudaMemcpyHostToDevice, *stream));
	checkCudaErrors(cudaMemcpyAsync(
			part_pos->z + idx_start_gpu, particles->z + idx_start_cpu,
			num_particles * sizeof(FPpart), cudaMemcpyHostToDevice, *stream));
	checkCudaErrors(cudaMemcpyAsync(
			part_pos->u + idx_start_gpu, particles->u + idx_start_cpu,
			num_particles * sizeof(FPpart), cudaMemcpyHostToDevice, *stream));
	checkCudaErrors(cudaMemcpyAsync(
			part_pos->v + idx_start_gpu, particles->v + idx_start_cpu,
			num_particles * sizeof(FPpart), cudaMemcpyHostToDevice, *stream));
	checkCudaErrors(cudaMemcpyAsync(
			part_pos->w + idx_start_gpu, particles->w + idx_start_cpu,
			num_particles * sizeof(FPpart), cudaMemcpyHostToDevice, *stream));
	checkCudaErrors(cudaMemcpyAsync(
			part_pos->q + idx_start_gpu, particles->q + idx_start_cpu,
			num_particles * sizeof(FPinterp), cudaMemcpyHostToDevice, *stream));
}

/**
 * Copy back num_particles number of particles from gpu memory to cpu
 */
void particles_positions_copy_to_host(
	cudaStream_t* stream,
	struct particles_positions_gpu* part_pos,
	struct particles* particles,
	size_t idx_start_cpu,
	size_t idx_start_gpu, 
	size_t num_particles
	) {

	checkCudaErrors(cudaMemcpyAsync(
			particles->x + idx_start_cpu, part_pos->x + idx_start_gpu,
			num_particles * sizeof(FPpart), cudaMemcpyDeviceToHost, *stream));
	checkCudaErrors(cudaMemcpyAsync(
			particles->y + idx_start_cpu, part_pos->y + idx_start_gpu,
			num_particles * sizeof(FPpart), cudaMemcpyDeviceToHost, *stream));
	checkCudaErrors(cudaMemcpyAsync(
			particles->z + idx_start_cpu, part_pos->z + idx_start_gpu,
			num_particles * sizeof(FPpart), cudaMemcpyDeviceToHost, *stream));
	checkCudaErrors(cudaMemcpyAsync(
			particles->u + idx_start_cpu, part_pos->u + idx_start_gpu,
			num_particles * sizeof(FPpart), cudaMemcpyDeviceToHost, *stream));
	checkCudaErrors(cudaMemcpyAsync(
			particles->v + idx_start_cpu, part_pos->v + idx_start_gpu,
			num_particles * sizeof(FPpart), cudaMemcpyDeviceToHost, *stream));
	checkCudaErrors(cudaMemcpyAsync(
			particles->w + idx_start_cpu, part_pos->w + idx_start_gpu,
			num_particles * sizeof(FPpart), cudaMemcpyDeviceToHost, *stream));
	checkCudaErrors(cudaMemcpyAsync(
			particles->q + idx_start_cpu, part_pos->q + idx_start_gpu,
			num_particles * sizeof(FPinterp), cudaMemcpyDeviceToHost, *stream));
}

/**
 *  Allocate and copy particle info to the device.
 */
void particles_info_alloc_and_copy_to_device(
		struct particles_info_gpu* part_info_gpu,
		struct particles_info_gpu** part_info_gpu_ptr,
		struct particles* particles_cpu) {
	// Copy data from the host object to the cpu object, the objects look
	// different, hence the manual copy.
	part_info_gpu->species_ID = particles_cpu->species_ID;
	part_info_gpu->npmax = particles_cpu->npmax;
	part_info_gpu->nop = particles_cpu->nop;
	part_info_gpu->NiterMover = particles_cpu->NiterMover;
	part_info_gpu->n_sub_cycles = particles_cpu->n_sub_cycles;
	part_info_gpu->npcel = particles_cpu->npcel;
	part_info_gpu->npcelx = particles_cpu->npcelx;
	part_info_gpu->npcely = particles_cpu->npcely;
	part_info_gpu->npcelz = particles_cpu->npcelz;
	part_info_gpu->qom = particles_cpu->qom;
	part_info_gpu->u0 = particles_cpu->u0;
	part_info_gpu->v0 = particles_cpu->v0;
	part_info_gpu->w0 = particles_cpu->w0;
	part_info_gpu->uth = particles_cpu->uth;
	part_info_gpu->vth = particles_cpu->vth;
	part_info_gpu->wth = particles_cpu->wth;

	checkCudaErrors(cudaMalloc(part_info_gpu_ptr, sizeof(particles_info_gpu)));
	checkCudaErrors(cudaMemcpy(*part_info_gpu_ptr, part_info_gpu,
														 sizeof(particles_info_gpu),
														 cudaMemcpyHostToDevice));
}

void particles_info_dealloc_device(
		struct particles_info_gpu* part_info_gpu_ptr) {
	checkCudaErrors(cudaFree(part_info_gpu_ptr));
}

/**
 * Update all particles of a species on gpu in batches of size batchsize.
 * Batchsize should be a multiple of 512
 *
 */
int batch_update_particles(cudaStream_t* stream, 
							struct particles* part_cpu,
							struct particles_positions_gpu* part_gpu,
							struct particles_positions_gpu* part_gpu_ptr,
							struct particles_info_gpu* part_info_gpu_ptr,
							struct EMfield* field_gpu_ptr,
							struct grid* grd_gpu_ptr,
							struct interpDensSpecies* ids_gpu_ptr,
							struct parameters* param_cpu,
							struct parameters* param_gpu_ptr, int batchsize) 
	{
	// TODO make entire function run on a cuda stream

	int offset_cpu = 0;
	int offset_gpu = 0;
	int current_batch;

	int batches = 0;

	while (offset_cpu < part_cpu->nop) {
		if (part_cpu->nop - offset_cpu >= batchsize) {
			current_batch = batchsize;
		} else { 
			current_batch = part_cpu->nop - offset_cpu;
		}

		particles_positions_copy_to_device(
			stream, 
			part_cpu, 
			part_gpu, 
			offset_cpu,
			offset_gpu,
			current_batch
			);
		move_and_interpolate<<<
			(current_batch + param_cpu->threads_per_block - 1) / param_cpu->threads_per_block,
			param_cpu->threads_per_block, 
			0, 
			*stream
			>>>(
			part_gpu_ptr, 
			part_info_gpu_ptr, 
			field_gpu_ptr, 
			grd_gpu_ptr,
			param_gpu_ptr, 
			ids_gpu_ptr, 
			current_batch
			);
		checkCudaErrors(cudaPeekAtLastError());

		particles_positions_copy_to_host(
			stream, 
			part_gpu, 
			part_cpu, 
			offset_cpu,
			offset_gpu,
			current_batch
			);

		offset_cpu += current_batch;
		batches++;
	}

	return batches;
}



__global__ void move_and_interpolate(struct particles_positions_gpu* part_pos,
																		 particles_info_gpu* part_info,
																		 struct EMfield* field, struct grid* grd,
																		 struct parameters* param,
																		 struct interpDensSpecies* ids_gpu,
																		 int num_particles) {
	int local_index = blockIdx.x * blockDim.x + threadIdx.x;

	if (local_index >= num_particles) {
		return;
	}

	mover_PC_gpu(part_pos, part_info, field, grd, param);
	interpP2G_gpu(part_pos, part_info, ids_gpu, grd);
}

__device__ void mover_PC_gpu(struct particles_positions_gpu* part_pos,
														 particles_info_gpu* part_info,
														 struct EMfield* field, struct grid* grd,
														 struct parameters* param) {
	int local_index = blockIdx.x * blockDim.x + threadIdx.x;

	// auxiliary variables
	FPpart dt_sub_cycling = (FPpart)param->dt / ((float)part_info->n_sub_cycles);
	FPpart dto2 = .5 * dt_sub_cycling, qomdt2 = part_info->qom * dto2 / param->c;
	FPpart omdtsq, denom, ut, vt, wt, udotb;

	// local (to the particle) electric and magnetic field
	FPfield Exl = 0.0, Eyl = 0.0, Ezl = 0.0, Bxl = 0.0, Byl = 0.0, Bzl = 0.0;

	// interpolation densities
	int ix, iy, iz;
	FPfield weight[2][2][2];
	FPfield xi[2], eta[2], zeta[2];

	// intermediate particle position and velocity
	FPpart xptilde, yptilde, zptilde, uptilde, vptilde, wptilde;

	// start subcycling
	for (int i_sub = 0; i_sub < part_info->n_sub_cycles; i_sub++) {
		// move each particle with new fields
		xptilde = part_pos->x[local_index];
		yptilde = part_pos->y[local_index];
		zptilde = part_pos->z[local_index];
		// calculate the average velocity iteratively
		for (int innter = 0; innter < part_info->NiterMover; innter++) {
			// interpolation G-->P
			ix = 2 + int((part_pos->x[local_index] - grd->xStart) * grd->invdx);
			iy = 2 + int((part_pos->y[local_index] - grd->yStart) * grd->invdy);
			iz = 2 + int((part_pos->z[local_index] - grd->zStart) * grd->invdz);

			// calculate weights
			xi[0] = part_pos->x[local_index] -
							grd->XN_flat[get_idx(ix - 1, iy, iz, grd->nyn, grd->nzn)];
			eta[0] = part_pos->y[local_index] -
							 grd->YN_flat[get_idx(ix, iy - 1, iz, grd->nyn, grd->nzn)];
			zeta[0] = part_pos->z[local_index] -
								grd->ZN_flat[get_idx(ix, iy, iz - 1, grd->nyn, grd->nzn)];
			xi[1] = grd->XN_flat[get_idx(ix, iy, iz, grd->nyn, grd->nzn)] -
							part_pos->x[local_index];
			eta[1] = grd->YN_flat[get_idx(ix, iy, iz, grd->nyn, grd->nzn)] -
							 part_pos->y[local_index];
			zeta[1] = grd->ZN_flat[get_idx(ix, iy, iz, grd->nyn, grd->nzn)] -
								part_pos->z[local_index];

			//#pragma unroll
			for (int ii = 0; ii < 2; ii++)
				//#pragma unroll
				for (int jj = 0; jj < 2; jj++)
					//#pragma unroll
					for (int kk = 0; kk < 2; kk++)
						weight[ii][jj][kk] = xi[ii] * eta[jj] * zeta[kk] * grd->invVOL;

			// set to zero local electric and magnetic field
			Exl = 0.0, Eyl = 0.0, Ezl = 0.0, Bxl = 0.0, Byl = 0.0, Bzl = 0.0;

			for (int ii = 0; ii < 2; ii++)
				for (int jj = 0; jj < 2; jj++)
					for (int kk = 0; kk < 2; kk++) {
						Exl += weight[ii][jj][kk] *
									 field->Ex_flat[get_idx(ix - ii, iy - jj, iz - kk, grd->nyn,
																					grd->nzn)];
						Eyl += weight[ii][jj][kk] *
									 field->Ey_flat[get_idx(ix - ii, iy - jj, iz - kk, grd->nyn,
																					grd->nzn)];
						Ezl += weight[ii][jj][kk] *
									 field->Ez_flat[get_idx(ix - ii, iy - jj, iz - kk, grd->nyn,
																					grd->nzn)];
						Bxl += weight[ii][jj][kk] *
									 field->Bxn_flat[get_idx(ix - ii, iy - jj, iz - kk, grd->nyn,
																					 grd->nzn)];
						Byl += weight[ii][jj][kk] *
									 field->Byn_flat[get_idx(ix - ii, iy - jj, iz - kk, grd->nyn,
																					 grd->nzn)];
						Bzl += weight[ii][jj][kk] *
									 field->Bzn_flat[get_idx(ix - ii, iy - jj, iz - kk, grd->nyn,
																					 grd->nzn)];
					}

			// end interpolation
			omdtsq = qomdt2 * qomdt2 * (Bxl * Bxl + Byl * Byl + Bzl * Bzl);
			denom = 1.0 / (1.0 + omdtsq);
			// solve the position equation
			ut = part_pos->u[local_index] + qomdt2 * Exl;
			vt = part_pos->v[local_index] + qomdt2 * Eyl;
			wt = part_pos->w[local_index] + qomdt2 * Ezl;
			udotb = ut * Bxl + vt * Byl + wt * Bzl;
			// solve the velocity equation
			uptilde =
					(ut + qomdt2 * (vt * Bzl - wt * Byl + qomdt2 * udotb * Bxl)) * denom;
			vptilde =
					(vt + qomdt2 * (wt * Bxl - ut * Bzl + qomdt2 * udotb * Byl)) * denom;
			wptilde =
					(wt + qomdt2 * (ut * Byl - vt * Bxl + qomdt2 * udotb * Bzl)) * denom;
			// update position
			part_pos->x[local_index] = xptilde + uptilde * dto2;
			part_pos->y[local_index] = yptilde + vptilde * dto2;
			part_pos->z[local_index] = zptilde + wptilde * dto2;

		}  // end of iteration
		// update the final position and velocity
		part_pos->u[local_index] = 2.0 * uptilde - part_pos->u[local_index];
		part_pos->v[local_index] = 2.0 * vptilde - part_pos->v[local_index];
		part_pos->w[local_index] = 2.0 * wptilde - part_pos->w[local_index];
		part_pos->x[local_index] = xptilde + uptilde * dt_sub_cycling;
		part_pos->y[local_index] = yptilde + vptilde * dt_sub_cycling;
		part_pos->z[local_index] = zptilde + wptilde * dt_sub_cycling;

		//////////
		//////////
		////////// BC

		bc_particles(part_pos->x, part_pos->u, param->PERIODICX, grd->Lx,
								 local_index);
		bc_particles(part_pos->y, part_pos->v, param->PERIODICY, grd->Ly,
								 local_index);
		bc_particles(part_pos->z, part_pos->w, param->PERIODICZ, grd->Lz,
								 local_index);

	}  // end of subcycling
}

__device__ void bc_particles(FPpart* pos, FPpart* vel, bool periodic,
														 double len, int index) {
	if (pos[index] > len) {
		if (periodic) {  // PERIODIC
			pos[index] = pos[index] - len;
		} else {  // REFLECTING BC
			vel[index] = -vel[index];
			pos[index] = 2 * len - pos[index];
		}
	} else if (pos[index] < 0) {
		if (periodic) {  // PERIODIC
			pos[index] = pos[index] + len;
		} else {  // REFLECTING BC
			vel[index] = -vel[index];
			pos[index] = -pos[index];
		}
	}
}

/** Interpolation Particle --> Grid: This is for species */
__device__ void interpP2G_gpu(struct particles_positions_gpu* part_pos,
															particles_info_gpu* part_info,
															struct interpDensSpecies* ids_gpu,
															struct grid* grd) {
	// Index and helper variables
	int local_index = blockIdx.x * blockDim.x + threadIdx.x;

	// arrays needed for interpolation
	FPpart weight[2][2][2];
	FPpart xi[2], eta[2], zeta[2];

	// index of the cell
	int ix, iy, iz;

	// determine cell: can we change to int()? is it faster?
	ix = 2 + int(floor((part_pos->x[local_index] - grd->xStart) * grd->invdx));
	iy = 2 + int(floor((part_pos->y[local_index] - grd->yStart) * grd->invdy));
	iz = 2 + int(floor((part_pos->z[local_index] - grd->zStart) * grd->invdz));

	// distances from node
	xi[0] = part_pos->x[local_index] -
					grd->XN_flat[get_idx(ix - 1, iy, iz, grd->nyn, grd->nzn)];
	eta[0] = part_pos->y[local_index] -
					 grd->YN_flat[get_idx(ix, iy - 1, iz, grd->nyn, grd->nzn)];
	zeta[0] = part_pos->z[local_index] -
						grd->ZN_flat[get_idx(ix, iy, iz - 1, grd->nyn, grd->nzn)];
	xi[1] = grd->XN_flat[get_idx(ix, iy, iz, grd->nyn, grd->nzn)] -
					part_pos->x[local_index];
	eta[1] = grd->YN_flat[get_idx(ix, iy, iz, grd->nyn, grd->nzn)] -
					 part_pos->y[local_index];
	zeta[1] = grd->ZN_flat[get_idx(ix, iy, iz, grd->nyn, grd->nzn)] -
						part_pos->z[local_index];

	/**
	 * From here on we are modifying shared data with other threads.
	 * Shared memory might speed this up.
	 */

	// calculate the weights for different nodes
	for (int ii = 0; ii < 2; ii++)
		for (int jj = 0; jj < 2; jj++)
			for (int kk = 0; kk < 2; kk++)
				weight[ii][jj][kk] = part_pos->q[local_index] * xi[ii] * eta[jj] *
														 zeta[kk] * grd->invVOL;

	//////////////////////////
	// add charge density
	for (int ii = 0; ii < 2; ii++)
		for (int jj = 0; jj < 2; jj++)
			for (int kk = 0; kk < 2; kk++)
				atomicAdd(&ids_gpu->rhon_flat[get_idx((ix - ii), (iy - jj), (iz - kk),
																							grd->nyn, grd->nzn)],
									weight[ii][jj][kk] * grd->invVOL);

	atomic_add_pressure(part_pos->u[local_index], weight, ids_gpu->Jx_flat,
											grd->invVOL, ix, iy, iz, grd->nyn, grd->nzn);
	atomic_add_pressure(part_pos->v[local_index], weight, ids_gpu->Jy_flat,
											grd->invVOL, ix, iy, iz, grd->nyn, grd->nzn);
	atomic_add_pressure(part_pos->w[local_index], weight, ids_gpu->Jz_flat,
											grd->invVOL, ix, iy, iz, grd->nyn, grd->nzn);

	atomic_add_pressure(part_pos->u[local_index] * part_pos->u[local_index],
											weight, ids_gpu->pxx_flat, grd->invVOL, ix, iy, iz,
											grd->nyn, grd->nzn);
	atomic_add_pressure(part_pos->v[local_index] * part_pos->v[local_index],
											weight, ids_gpu->pyy_flat, grd->invVOL, ix, iy, iz,
											grd->nyn, grd->nzn);
	atomic_add_pressure(part_pos->w[local_index] * part_pos->w[local_index],
											weight, ids_gpu->pzz_flat, grd->invVOL, ix, iy, iz,
											grd->nyn, grd->nzn);

	atomic_add_pressure(part_pos->u[local_index] * part_pos->v[local_index],
											weight, ids_gpu->pxy_flat, grd->invVOL, ix, iy, iz,
											grd->nyn, grd->nzn);
	atomic_add_pressure(part_pos->u[local_index] * part_pos->w[local_index],
											weight, ids_gpu->pxz_flat, grd->invVOL, ix, iy, iz,
											grd->nyn, grd->nzn);
	atomic_add_pressure(part_pos->v[local_index] * part_pos->w[local_index],
											weight, ids_gpu->pyz_flat, grd->invVOL, ix, iy, iz,
											grd->nyn, grd->nzn);
}

__device__ void atomic_add_pressure(FPpart part_dat, FPpart weight[2][2][2],
																		FPinterp* ids_arr, FPfield invVOL, int ix,
																		int iy, int iz, int stride_y,
																		int stride_z) {
	FPpart temp[2][2][2];
	// #pragma unroll
	for (int ii = 0; ii < 2; ii++)
		//#pragma unroll
		for (int jj = 0; jj < 2; jj++)
			//#pragma unroll
			for (int kk = 0; kk < 2; kk++)
				temp[ii][jj][kk] = part_dat * weight[ii][jj][kk];

	//#pragma unroll
	for (int ii = 0; ii < 2; ii++)
		//#pragma unroll
		for (int jj = 0; jj < 2; jj++)
			//#pragma unroll
			for (int kk = 0; kk < 2; kk++)
				atomicAdd(&ids_arr[get_idx((ix - ii), (iy - jj), (iz - kk), stride_y,
																	 stride_z)],
									temp[ii][jj][kk] * invVOL);
}
