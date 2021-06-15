/** A mixed-precision implicit Particle-in-Cell simulator for heterogeneous systems **/

// Allocator for 2D, 3D and 4D array: chain of pointers
#include "Alloc.h"

// Precision: fix precision for different quantities
#include "PrecisionTypes.h"
// Simulation Parameter - structure
#include "Parameters.h"
// Grid structure
#include "Grid.h"
// Interpolated Quantities Structures
#include "InterpDensSpecies.h"
#include "InterpDensNet.h"
#include "InterpDens_aux.h"

// Field structure
#include "EMfield.h" // Just E and Bn
//#include "EMfield_aux.h" // Bc, Phi, Eth, D

// Particles structure
#include "Particles.h"
//#include "Particles_aux.h" // Needed only if dointerpolation on GPU - avoid reduction on GPU

// solvers
#include "MaxwellSolver.h"

// mover
#include "Mover.h"

// Initial Condition
#include "IC.h"
// Boundary Conditions
#include "BC.h"
// Smoothing
#include "Smoothing.h"
// timing
#include "Timing.h"
// Read and output operations
#include "RW_IO.h"

#include <omp.h>
#include <mpi.h>
#include <sys/stat.h>
#include <stdexcept>

#include "mpi_comm.h"

#if defined(USE_CATALYST) && 0
#include "Adaptor.h"
#endif


// ====================================================== //
// Local function declarations

double timer(
    double *mean, 
    double *variance, 
    double *cycle,
    double start_time,
    long count);


// ====================================================== //
// Main


int main(int argc, char **argv){


    // ====================================================== //
    // Init MPI 
	int mpi_thread_support;
	int mpi_rank, mpi_size;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_thread_support);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    std::cout << "Total number of cores: " << omp_get_max_threads() << std::endl;
    

    // ====================================================== //
    // Read the inputfile and fill the param structure
    // Read the input file name from command line

    parameters param;
    readInputFile(&param,argc,argv);
    if(!mpi_rank){
        printParameters(&param);
        saveParameters(&param);
        struct stat sb;
        if (stat(param.RestartDirName.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
            std::cout << "Output directory " + param.RestartDirName + " exists." << std::endl;
        }
        else {
            std::cerr << "Output directory " + param.RestartDirName + " does not exists." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // ====================================================== //
    // Declare variables and alloc memory

        // Set-up the grid information
    grid grd;
    setGrid(&param, &grd);
    
    // Allocate Fields
    EMfield field;
    field_allocate(&grd,&field);
    EMfield_aux field_aux;
    field_aux_allocate(&grd,&field_aux);
    
    // Allocate Interpolated Quantities
    // per species
    interpDensSpecies *ids = new interpDensSpecies[param.ns];
    for (int is=0; is < param.ns; is++)
        interp_dens_species_allocate(&grd,&ids[is],is);
    // Net densities
    interpDensNet idn;
    interp_dens_net_allocate(&grd,&idn);
    // Hat densities
    interpDens_aux id_aux;
    interp_dens_aux_allocate(&grd,&id_aux);
    
    
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

#if defined(USE_CATALYST) && 0
    if (!mpi_rank)
        printf("init catalyst\n");
        Adaptor::Initialize("../scripts/image.py",
      	                    0,
      	                    0,
      	                    0,
      	                    grd.nxn,
      	                    grd.nyn,
      	                    grd.nzn,
      	                    grd.dx,
      	                    grd.dy,
      	                    grd.dz);
#endif

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
        particle_allocate(&param,&part[is],is);

    mpi_broadcast_field(&grd, &field);
    for(int is=0; is<param.ns; is++){
        mpi_scatter_particles(&part_global[is], &part[is]);
    }

    // Dealloc global particles array, it is no longer needed
    if(!mpi_rank){
        for (int is=0; is < param.ns; is++)
            particle_deallocate(&part_global[is]);
    }



    // ====================================================== //
    // Timing variables
    double iStart = MPI_Wtime();
    double time0 = iStart;
    // avg, variance and last cycle exec time for mover, interpolation, field solver and io, respectively
    double average[5] = {0.,0.,0.,0.,0.};
    double variance[5] = {0.,0.,0.,0.,0.};
    double cycle_time[5] = {0.,0.,0.,0.,0.};
    
    // **********************************************************//
    // **** Start the Simulation!  Cycle index start from 1  *** //
    // **********************************************************//
    for (int cycle = param.first_cycle_n; cycle < (param.first_cycle_n + param.ncycles); cycle++) {

        // ====================================================== //
        // implicit mover

        if(!mpi_rank){
            std::cout << std::endl;
            std::cout << "***********************" << std::endl;
            std::cout << "   cycle = " << cycle << std::endl;
            std::cout << "***********************" << std::endl;

            std::cout << "*** Particle Mover ***" << std::endl;
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

        // #pragma omp parallel for // only if use mover_PC_V
        for (int is=0; is < param.ns; is++){
            mover_PC(&part[is],&field,&grd,&param);
            //mover_PC_V(&part[is],&field,&grd,&param);
            //mover_interp(&part[is], &field, &ids[is],&grd, &param);

        }

        time0 = timer(&average[1], &variance[1], &cycle_time[1], time0, cycle);

        // ====================================================== //
        // interpolation particle to grid
        if(!mpi_rank)
            std::cout << "*** INTERPOLATION P2G ***" << std::endl;

        // set to zero the densities - needed for interpolation
        setZeroDensities(&idn,ids,&grd,param.ns);

        // interpolate species: MAXIMUM parallelism is number of species
        #pragma omp parallel for
        for (int is=0; is < param.ns; is++){
            interpP2G(&part[is],&ids[is],&grd);

           // apply BC to interpolated densities
            applyBCids(&ids[is],&grd,&param);
        }

        // sum over species
        sumOverSpecies(&idn,ids,&grd,param.ns);

        // MPI communicate densities to master, both net densities and for species
		mpi_reduce_densities(&grd, &idn);
		for(int is=0; is<param.ns; is++){
			mpi_reduce_densities(&grd, &ids[is]);
		}

        // Update timer for interpolation
        time0 = timer(&average[2], &variance[2], &cycle_time[2], time0, cycle);

        // ====================================================== //
        // From here, master calculates new EM field. slaves idle

        if(!mpi_rank){
            // interpolate charge density from center to node
            applyBCscalarDensN(idn.rhon,&grd,&param);
            interpN2Cinterp(idn.rhoc,idn.rhon, &grd);

            // ====================================================== //
            // Maxwell solver

            // calculate hat functions rhoh and Jxh, Jyh, Jzh
            calculateHatDensities(&id_aux, &idn, ids, &field, &grd, &param);

            //  Poisson correction
            if (param.PoissonCorrection)
                divergenceCleaning(&grd,&field_aux,&field,&idn,&param);
            
            calculateE(&grd,&field_aux,&field,&id_aux,ids,&param);
            calculateB(&grd,&field_aux,&field,&param);

        }

        // broadcast EM field data from master process to slaves
        mpi_broadcast_field(&grd, &field);
        // Update timer for field solver
        time0 = timer(&average[3], &variance[3], &cycle_time[3], time0, cycle);
            
        // ====================================================== //
        // IO
        if(!mpi_rank){
            // write E, B, rho to disk
            if (cycle%param.FieldOutputCycle==0){
                VTK_Write_Vectors_Binary(cycle, &grd,&field, &param);
                VTK_Write_Scalars_Binary(cycle, &grd,ids,&idn, &param);
            }
        }

        if (cycle % param.ParticlesOutputCycle == 0) {
             HDF5_Write_Particles(cycle, part, &param);
        }

#if defined(USE_CATALYST) && 0
        printf("CoProcess\n");
        Adaptor::CoProcess(param.dt*cycle, cycle, field.Bxn, field.Byn, field.Bzn, ids[0].rhon, ids[1].rhon);
#endif

        // Update timer for io
        time0 = timer(&average[4], &variance[4], &cycle_time[4], time0, cycle);

        if(!mpi_rank)
            std::cout << "Timing Cycle " << cycle << " : " << cycle_time[0] << " " << cycle_time[1] << " " << cycle_time[2] << " " << cycle_time[3] << " " << cycle_time[4] << std::endl;
    }  // end of one PIC cycle
    
    /// Release the resources
    // deallocate field
    grid_deallocate(&grd);
    field_deallocate(&grd,&field);
    // interp
    interp_dens_net_deallocate(&grd,&idn);
    interp_dens_aux_deallocate(&grd,&id_aux);
    
    // Deallocate interpolated densities and particles
    for (int is=0; is < param.ns; is++){
        interp_dens_species_deallocate(&grd,&ids[is]);
        particle_deallocate(&part[is]);
    }
    
    if(!mpi_rank){
        // stop timer
        double iElaps = MPI_Wtime() - iStart;

        // Print timing of simulation
        std::cout << std::endl;
        std::cout << "******************************************************" << std::endl;
        std::cout << "   Tot. Simulation Time (s) = " << iElaps << std::endl;
        std::cout << "   Part. Sort / Cycle   (s) = " << average[0] << "+-" << sqrt(variance[0] / (param.ncycles - 1)) << std::endl;
        std::cout << "   Mover Time / Cycle   (s) = " << average[1] << "+-" << sqrt(variance[1] / (param.ncycles - 1)) << std::endl;
        std::cout << "   Interp. Time / Cycle (s) = " << average[2] << "+-" << sqrt(variance[2] / (param.ncycles - 1)) << std::endl;
        std::cout << "   Field Time / Cycle   (s) = " << average[3] << "+-" << sqrt(variance[3] / (param.ncycles - 1)) << std::endl;
        std::cout << "   IO Time / Cycle      (s) = " << average[4] << "+-" << sqrt(variance[4] / (param.ncycles - 1)) << std::endl;
        std::cout << "******************************************************" << std::endl;
    }

#if defined(USE_CATALYST) && 0
    Adaptor::Finalize();
#endif

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
