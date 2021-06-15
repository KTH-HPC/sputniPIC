/** A mixed-precision implicit Particle-in-Cell simulator for heterogeneous
 * systems **/

#include <mpi.h>
#include <mpi_comm.h>

#include <omp.h>
#include <sys/stat.h>
#include <stdexcept>


// Allocator for 2D, 3D and 4D array: chain of pointers
#include "Alloc.h"

// Precision: fix precision for different quantities
#include "PrecisionTypes.h"
// Simulation Parameter - structure
#include "Parameters.h"
#include "gpu/Parameters_gpu.h"
// Grid structure
#include "Grid.h"
#include "gpu/Grid_gpu.h"
// Interpolated Quantities Structures
#include "gpu/InterpDensNet_gpu.cuh"
#include "gpu/InterpDensSpecies_gpu.h"

// Field structure
#include "EMfield.h"          // Just E and Bn
#include "gpu/EMfield_gpu.h"  // Just E and Bn

// Particles structure
#include "Particles.h"
#include "gpu/Particles_gpu.cuh"

// solvers
#include "MaxwellSolver.h"

// mover
#include "Mover.h"

// Initial Condition
#include "IC.h"
// Boundary Conditions
#include "BC.h"
#include "gpu/BC_gpu.cuh"
// Smoothing
#include "Smoothing.h"
// timing
#include "Timing.h"
// Read and output operations
#include "RW_IO.h"

// Cuda memcheck and particle batching helper
#include "gpu/cuda_helper.h"


double timer(
    double *mean, 
    double *variance, 
    double *cycle,
    double start_time,
    long count);

// ====================================================== //
// Main


int main(int argc, char** argv) {

	int mpi_thread_support;
	int mpi_rank, mpi_size;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_thread_support);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // ====================================================== //
    // Read the inputfile and fill the param structure
    // Read the input file name from command line

    parameters param;
    readInputFile(&param,argc,argv);
    if(!mpi_rank){
        printParameters(&param);
        saveParameters(&param);
    }

	if(!mpi_rank){
		struct stat sb;
		if (stat(param.RestartDirName.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
			std::cout << "Output directory " + param.RestartDirName + " exists." << std::endl;
		}
		else {
			throw std::runtime_error("Output directory " + param.RestartDirName + " does not exists.");
		}
	}

	// ====================================================== //
    // Declare variables and alloc memory

	// Set-up the grid information
	grid grd;
	setGrid(&param, &grd);

	// Allocate Fields
	EMfield field;
	field_allocate(&grd, &field);
	EMfield_aux field_aux;
	field_aux_allocate(&grd, &field_aux);

	// Allocate Interpolated Quantities
	// per species
	interpDensSpecies* ids = new interpDensSpecies[param.ns];
	for (int is = 0; is < param.ns; is++) {
		interp_dens_species_allocate(&grd, &ids[is], is);
	}

	// Net densities
	interpDensNet idn;
	interp_dens_net_allocate(&grd, &idn);

	// Hat densities
	interpDens_aux id_aux;
	interp_dens_aux_allocate(&grd, &id_aux);

    // Allocate Particles
    particles *part = new particles[param.ns];
    particles *part_global = new particles[param.ns];
    
    // allocation for global particles
    if(!mpi_rank){
        for (int is=0; is < param.ns; is++){
            particle_allocate(&param, &part_global[is], is);
        }
    }

    // ====================================================== //
    // Initialization global system

    if(!mpi_rank){
        initGEM(&param,&grd,&field,&field_aux,part_global,ids);
        //initUniform(&params_global,&grd,&field,&field_aux,part_global,ids);
    }

    // ====================================================== //
    // Distribute system to slave processors.
    // We do particle decomposition, shared domain. 

    // Set number of particles per species for local processes
	for (int i = 0; i < param.ns; i++){
		//number of particles localy is global number/num mpi processes
		param.np[i] /= mpi_size;
		// Maximum number of particles is also divided by num mpi processes,
		// since we do particle decomposition
		param.npMax[i] /= mpi_size;
	}

    // allocation of local particles
    for (int is=0; is < param.ns; is++)
		// Allocate local particles in pinned memory
        particle_allocate_host(&param,&part[is],is);

    mpi_broadcast_field(&grd, &field);
    for(int is=0; is<param.ns; is++){
        mpi_scatter_particles(&part_global[is], &part[is]);
    }

    // Dealloc global particles array, it is no longer needed
    if(!mpi_rank){
        for (int is=0; is < param.ns; is++)
            particle_deallocate(&part_global[is]);
    }


	// ********************** GPU ALLOCATION ********************//
	// Create the parameters on the GPU and copy the data to the device.

	// Get number of gpu devices
	int num_devices;
	checkCudaErrors(cudaGetDeviceCount(&num_devices));

	parameters* param_gpu_ptr[num_devices] = {nullptr};
	for (int device_id = 0; device_id < num_devices; device_id++) {
		checkCudaErrors(cudaSetDevice(device_id));
		param_alloc_and_copy_to_device(&param, &param_gpu_ptr[device_id]);
	}

	// Create the grid on the GPU and copy the data to the device.
	grid grid_gpu[num_devices];
	grid* grid_gpu_ptr[num_devices] = {nullptr};
	for (int device_id = 0; device_id < num_devices; device_id++) {
		checkCudaErrors(cudaSetDevice(device_id));
		grid_alloc_and_copy_to_device(&grd, &grid_gpu[device_id],
										&grid_gpu_ptr[device_id]);
	}

	// Create the electromagnetic field on the GPU and copy the data to the
	// device.
	EMfield field_gpu[num_devices];
	EMfield* field_gpu_ptr[num_devices] = {nullptr};
	for (int device_id = 0; device_id < num_devices; device_id++) {
		checkCudaErrors(cudaSetDevice(device_id));
		field_alloc_and_copy_to_device(&grd, 
										&field_gpu[device_id], 
										&field,
										&field_gpu_ptr[device_id]);
	}

	// Create interpolated quantities on the GPU and and copy the data to the
	// device.
	interpDensSpecies* ids_gpu = new interpDensSpecies[param.ns];        // array
	interpDensSpecies** ids_gpu_ptr = new interpDensSpecies*[param.ns];  // array
	for (int is = 0; is < param.ns; is++) {
		int device_id = (is + num_devices) % num_devices;
		checkCudaErrors(cudaSetDevice(device_id));
		// Init and copy data here.
		interp_DS_alloc_and_copy_to_device(&ids[is], &ids_gpu[is], &ids_gpu_ptr[is],
											&grd);
	}

	// Create particle information on the GPU and copy the data to the device.
	particles_info_gpu* part_info_gpu =
			new particles_info_gpu[param.ns];  // array
	particles_info_gpu** part_info_gpu_ptr =
			new particles_info_gpu*[param.ns];  // array

	particles_positions_gpu* part_positions_gpu =
			new particles_positions_gpu[param.ns];
	particles_positions_gpu** part_positions_gpu_ptr =
			new particles_positions_gpu*[param.ns];

	size_t np = 0;
	for (int is = 0; is < param.ns; is++) {
		np += part[is].nop;
	}

	size_t batch_size = get_appropriate_batch_size(&param);
	if (batch_size <= 0) {
		return -1;
	}
//  batch_size /= num_devices; /* size computed based on free mem of one device */

	// allocate field and copy to GPUs
	for (int is = 0; is < param.ns; is++) {
		int device_id = (is + num_devices) % num_devices;
		checkCudaErrors(cudaSetDevice(device_id));
		particles_info_alloc_and_copy_to_device(
			&part_info_gpu[is],
			&part_info_gpu_ptr[is], 
			&part[is]
			);
		particles_positions_alloc_device(
			&part_positions_gpu[is],
			&part_positions_gpu_ptr[is], 
			batch_size
			);
	}

	// create cuda streams for asynchronous workloading
	cudaStream_t* streams = new cudaStream_t[param.ns];
	for (int is = 0; is < param.ns; is++) {
		int device_id = (is + num_devices) % num_devices;
		checkCudaErrors(cudaSetDevice(device_id));
		checkCudaErrors(
				cudaStreamCreateWithFlags(&streams[is], cudaStreamNonBlocking));
	}

	if(mpi_rank == 0){

		#ifdef MEMCHECK
		std::cout << "CUDA return check: Enabled" << std::endl;
		#else
		std::cout << "CUDA return check: Disabled" << std::endl;
		#endif
		std::cout << "Total number of MPI ranks: " << mpi_size << std::endl;
		std::cout << "Number of cores per rank: " << omp_get_max_threads() << std::endl;
		std::cout << "Number of GPUs per rank: " << num_devices << std::endl;
		std::cout << "Threads Per Block of GPUs: " << param.threads_per_block << std::endl;
		std::cout << "Total number of particles: " << np*mpi_size << std::endl;
		std::cout << "Number of particles per rank: " << np << "; "
				<< (np * sizeof(FPpart) * 6) / (1024 * 1024) << " MB of data"
				<< std::endl;
		std::cout << "Allocating "
				<< (batch_size * param.ns * sizeof(FPpart) * 6) / (1024 * 1024)
				<< " MB of memory for particles on gpu" << std::endl;
		std::cout << "batch_size per species of " << batch_size << " ("
				<< (batch_size * sizeof(FPpart) * 6) / (1024 * 1024) << " MB)"
				<< std::endl;
	}


	// ====================================================== //
    // Timing variables
    double iStart = MPI_Wtime();
    double time0 = iStart;
    // avg, variance and last cycle exec time for sort, mover, interpolation, field solver and io, respectively
    double average[5] = {0.,0.,0.,0.,0.};
    double variance[5] = {0.,0.,0.,0.,0.};
    double cycle_time[5] = {0.,0.,0.,0.,0.};

	// **********************************************************//
	// **** Start the Simulation!  Cycle index start from 1  *** //
	// **********************************************************//
	for (int cycle = param.first_cycle_n;
			 cycle < (param.first_cycle_n + param.ncycles); cycle++) {


		if(mpi_rank == 0){
			std::cout << std::endl;
			std::cout << "***********************" << std::endl;
			std::cout << "   cycle = " << cycle << std::endl;
			std::cout << "***********************" << std::endl;

		}

		time0 = MPI_Wtime();

                if (param.sort) {
                  if (cycle % param.sort_every_n == 0 && cycle != 0) {
                    std::cout << "Sorting particles." << std::endl;
                    for (int is = 0; is < param.ns; is++) {
                      particle_sort(&param, &part[is], &grd);
                    }
                  }
                }

                // timer for sorting
		time0 = timer(&average[0], &variance[0], &cycle_time[0], time0, cycle);

		// set to zero the densities - needed for interpolation
		setZeroNetDensities(&idn, &grd);

		// blocking copy updated field from prev cycle
		for (int device_id = 0; device_id < num_devices; device_id++) {
			checkCudaErrors(cudaSetDevice(device_id));
			copy_from_host_to_device(&field_gpu[device_id], &field,
									grd.nxn * grd.nyn * grd.nzn);
		}

		// ====================================================== //
        // implicit mover

		// async set zero, move and interp
		for (int is = 0; is < param.ns; is++) {
			int device_id = (is + num_devices) % num_devices;

			checkCudaErrors(cudaSetDevice(device_id));
			setZeroSpeciesDensities_gpu(&streams[is], &grd, grid_gpu_ptr[device_id],
										&ids_gpu[is], is);
			
			// Move particles and interpolate to grid, on gpu in batches
			int b = batch_update_particles(
					&streams[is], &part[is], &part_positions_gpu[is],
					part_positions_gpu_ptr[is], part_info_gpu_ptr[is],
					field_gpu_ptr[device_id], grid_gpu_ptr[device_id], ids_gpu_ptr[is],
					&param, param_gpu_ptr[device_id], batch_size);

			if(mpi_rank == 0)
				std::cout << "***  MOVER  ITERATIONS = " << part[is].NiterMover 
					<< " - Species " << part[is].species_ID << " ***" 
					<< " on gpu " << device_id << " - " << b << " batches " 
					<< std::endl;
		
		}
		if(mpi_rank == 0)
			std::cout << "***********************" << std::endl;
		// Apply BCs
		for (int is = 0; is < param.ns; is++) {
			/**
			 * synchronization needed, interpolation needs to finish before we can
			 * apply boundary conditions
			 */
			int device_id = (is + num_devices) % num_devices;
			checkCudaErrors(cudaSetDevice(device_id));
			cudaStreamSynchronize(streams[is]);
			applyBCids_gpu(&streams[is], ids_gpu_ptr[is], grid_gpu_ptr[device_id],
										 &grd, &param);
		}

		// copy back quantities from GPUs
		for (int is = 0; is < param.ns; is++) {
			int device_id = (is + num_devices) % num_devices;
			checkCudaErrors(cudaSetDevice(device_id));
			interp_DS_copy_to_host(&ids_gpu[is], &ids[is], &grd);
		}

		// sum over species
		sumOverSpecies(&idn, ids, &grd, param.ns);

		// MPI communicate densities to root
		mpi_reduce_densities(&grd, &idn);
		for(int is=0; is<param.ns; is++){
			mpi_reduce_densities(&grd, &ids[is]);
		}

		// Update timer for mover+interpolation
		time0 = timer(&average[1], &variance[1], &cycle_time[1], time0, cycle);

		if(mpi_rank == 0){

			// interpolate charge density from center to node
			applyBCscalarDensN(idn.rhon, &grd, &param);
			interpN2Cinterp(idn.rhoc, idn.rhon, &grd);
			// calculate hat functions rhoh and Jxh, Jyh, Jzh
			calculateHatDensities(&id_aux, &idn, ids, &field, &grd, &param);

            // ====================================================== //
            // Maxwell solver

			//  Poisson correction
			if (param.PoissonCorrection)
				divergenceCleaning(&grd, &field_aux, &field, &idn, &param);

			calculateE(&grd, &field_aux, &field, &id_aux, ids, &param);
			calculateB(&grd, &field_aux, &field, &param);

		}

        // broadcast EM field data from master process to slaves
        mpi_broadcast_field(&grd, &field);
        // Update timer for field solver
		time0 = timer(&average[3], &variance[3], &cycle_time[3], time0, cycle);


		// ====================================================== //
        // IO
		if (!mpi_rank && cycle % param.FieldOutputCycle == 0) {
			// write E, B, rho to disk
			VTK_Write_Vectors(cycle, &grd, &field, &param);
			VTK_Write_Scalars(cycle, &grd, ids, &idn, &param);
		}

                if (cycle % param.ParticlesOutputCycle == 0) {
                    HDF5_Write_Particles(cycle, part, &param);
                }

        // Update timer for io
        time0 = timer(&average[4], &variance[4], &cycle_time[4], time0, cycle);

        if(!mpi_rank)
            std::cout << "Timing Cycle " << cycle << " : " << cycle_time[0] << " " << cycle_time[1] << " " << cycle_time[2] << " " << cycle_time[3] << " " << cycle_time[4] << std::endl;

	}  // end of one PIC cycle

	/// Release the resources
	// deallocate field
	grid_deallocate(&grd);
	field_deallocate(&grd, &field);
	// interp
	interp_dens_net_deallocate(&grd, &idn);
	interp_dens_aux_deallocate(&grd, &id_aux);

	// Deallocate interpolated densities and particles
	for (int is = 0; is < param.ns; is++) {
		int device_id = (is + num_devices) % num_devices;
		checkCudaErrors(cudaSetDevice(device_id));

		interp_dens_species_deallocate(&grd, &ids[is]);
		particle_deallocate_host(&part[is]);
	}

	if(mpi_rank == 0)
		std::cout << "DEALLOCATED CPU RESOURCES" << std::endl;

	// Deallocate GPU resources
	// streams
	for (int is = 0; is < param.ns; is++) {
		checkCudaErrors(cudaStreamDestroy(streams[is]));
	}

	for (int device_id = 0; device_id < num_devices; device_id++) {
		// Parameters
		checkCudaErrors(cudaFree(param_gpu_ptr[device_id]));
		// Grid
		grid_deallocate_device(&grid_gpu[device_id], grid_gpu_ptr[device_id]);
		// EMfield
		field_deallocate_device(&field_gpu[device_id], field_gpu_ptr[device_id]);
	}

	// interpDensSpecies
	for (int is = 0; is < param.ns; is++) {
		int device_id = (is + num_devices) % num_devices;
		checkCudaErrors(cudaSetDevice(device_id));

		interp_DS_dealloc_device(&ids_gpu[is], ids_gpu_ptr[is]);
		particles_info_dealloc_device(part_info_gpu_ptr[is]);
		particles_positions_dealloc_device(&part_positions_gpu[is], part_positions_gpu_ptr[is]);
	}


	delete[] ids_gpu;
	delete[] ids_gpu_ptr;
	delete[] part_info_gpu;
	delete[] part_info_gpu_ptr;
	delete[] part_positions_gpu;
	delete[] part_positions_gpu_ptr;

	if(mpi_rank == 0)
		std::cout << "DEALLOCATED GPU RESOURCES" << std::endl;


    if(!mpi_rank){
        // stop timer
        double iElaps = MPI_Wtime() - iStart;

        // Print timing of simulation
        std::cout << std::endl;
        std::cout << "******************************************************" << std::endl;
        std::cout << "   Tot. Simulation Time (s) = " << iElaps << std::endl;
        std::cout << "   Paticle sort / Cycle (s) = " << average[0] << "+-" << sqrt(variance[0] / (param.ncycles - 1)) << std::endl;
        std::cout << "   Mover/Interp / Cycle (s) = " << average[1] << "+-" << sqrt(variance[1] / (param.ncycles - 1)) << std::endl;
        // std::cout << "   Interp. Time / Cycle (s) = " << average[1] << "+-" << sqrt(variance[1] / (param.ncycles - 1)) << std::endl;
        std::cout << "   Field Time / Cycle   (s) = " << average[3] << "+-" << sqrt(variance[3] / (param.ncycles - 1)) << std::endl;
        std::cout << "   IO Time / Cycle      (s) = " << average[4] << "+-" << sqrt(variance[4] / (param.ncycles - 1)) << std::endl;
        std::cout << "******************************************************" << std::endl;
    }

    MPI_Finalize();
    // exit
    return 0;
}

inline void update_statistics(
    double *mean, 
    double *variance, 
    double new_value, 
    long count
    )
{
    double delta = new_value - *mean;
    *mean += delta / (double)count;
    double delta2 = new_value - *mean;
    *variance += delta * delta2;
}

double timer(
    double *mean, 
    double *variance, 
    double *cycle,
    double start_time,
    long count)
{
    double t = MPI_Wtime();
    *cycle = t - start_time;
    update_statistics(mean, variance, *cycle, count);
    return t;
}

