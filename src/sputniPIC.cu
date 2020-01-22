/** A mixed-precision implicit Particle-in-Cell simulator for heterogeneous
 * systems **/

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
#include "EMfield.h"  // Just E and Bn
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

#include "gpu/cuda_helper.h"

int main(int argc, char **argv) {
  // Read the inputfile and fill the param structure
  parameters param;
  // Read the input file name from command line
  readInputFile(&param, argc, argv);
  printParameters(&param);
  saveParameters(&param);

  // Timing variables
  double iStart = cpuSecond();
  double iMover, iInterp, iField, iIO, eMover = 0.0, eInterp = 0.0,
                                       eField = 0.0, eIO = 0.0;
  double dMover, dInterp, dField;

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
  interpDensSpecies *ids = new interpDensSpecies[param.ns];
  for (int is = 0; is < param.ns; is++)
    interp_dens_species_allocate(&grd, &ids[is], is);
  // Net densities
  interpDensNet idn;
  interp_dens_net_allocate(&grd, &idn);
  // Hat densities
  interpDens_aux id_aux;
  interp_dens_aux_allocate(&grd, &id_aux);

  // Allocate Particles
  particles *part = new particles[param.ns];
  // allocation
  for (int is = 0; is < param.ns; is++) {
    particle_allocate_host(&param, &part[is], is);
  }

  // Initialization
  initGEM(&param, &grd, &field, &field_aux, part, ids);
  // initUniform(&param,&grd,&field,&field_aux,part,ids);
  // ********************** GPU ALLOCATION ********************//
  // Create the parameters on the GPU and copy the data to the device.
  parameters* param_gpu_ptr = nullptr;
  param_alloc_and_copy_to_device(&param, &param_gpu_ptr);

  // Create the grid on the GPU and copy the data to the device.
  grid grid_gpu;
  grid* grid_gpu_ptr = nullptr;
  grid_alloc_and_copy_to_device(&grd, &grid_gpu, &grid_gpu_ptr);
  
  // Create the electromagnetic field on the GPU and copy the data to the device.
  EMfield field_gpu;
  EMfield* field_gpu_ptr = nullptr;
  field_alloc_and_copy_to_device(&grd, &field_gpu, &field, &field_gpu_ptr);
  
  // Create interpolated quantities on the GPU and and copy the data to the device.
  interpDensSpecies* ids_gpu = new interpDensSpecies[param.ns]; //array
  interpDensSpecies** ids_gpu_ptr = new interpDensSpecies*[param.ns]; //array
  for (int is = 0; is < param.ns; is++)
  {
      // Init and copy data here.
      interp_DS_alloc_and_copy_to_device(&ids[is], &ids_gpu[is], &ids_gpu_ptr[is], &grd);
  }

  interpDensNet idn_gpu;
  interpDensNet* idn_gpu_ptr = nullptr;
  interp_DNet_alloc_and_copy_to_device(&idn, &idn_gpu, &idn_gpu_ptr, &grd);

  // Create particle information on the GPU and copy the data to the device.
  particles_info_gpu* part_info_gpu = new particles_info_gpu[param.ns]; // array
  particles_info_gpu** part_info_gpu_ptr = new particles_info_gpu*[param.ns]; // array

  particles_positions_gpu* part_positions_gpu = new particles_positions_gpu[param.ns];
  particles_positions_gpu** part_positions_gpu_ptr = new particles_positions_gpu*[param.ns];

  int np = 0;
  for(int is=0; is<param.ns; is++){
      np += part[is].nop;
  }
  int batch_size = get_appropriate_batch_size(param.ns);
  if(batch_size <= 0){
      return -1;
  }

  std::cout << "Total number of particles: " << np << "; " << (np*sizeof(FPpart)*6)/(1024*1024) << " MB of data."  << std::endl;
  std::cout << "Allocating " << (batch_size*param.ns*sizeof(FPpart)*6)/(1024*1024) << " MB of memory for particles on gpu" << std::endl;
  std::cout << "(batch_size of " << (batch_size*sizeof(FPpart)*6)/(1024*1024) << " MB)" << std::endl;
  for (int is = 0; is < param.ns; is++)
  {
      particles_info_alloc_and_copy_to_device(&part_info_gpu[is], &part_info_gpu_ptr[is], &part[is]);
      particles_positions_alloc_device(&part_positions_gpu[is], &part_positions_gpu_ptr[is], batch_size);
  }  

  // create cuda streams for asynchronous workloading 
  cudaStream_t* streams = new cudaStream_t[param.ns];
  for(int is=0; is<param.ns; is++){

      checkCudaErrors(cudaStreamCreateWithFlags(&streams[is], cudaStreamNonBlocking));
  }



  // **********************************************************//
  // **** Start the Simulation!  Cycle index start from 1  *** //
  // **********************************************************//
  for (int cycle = param.first_cycle_n;
       cycle < (param.first_cycle_n + param.ncycles); cycle++) {
    std::cout << std::endl;
    std::cout << "***********************" << std::endl;
    std::cout << "   cycle = " << cycle << std::endl;
    std::cout << "***********************" << std::endl;

    // set to zero the densities - needed for interpolation
//    setZeroDensities(&idn, ids, &grd, param.ns);
    copy_from_host_to_device(&field_gpu, &field, grd.nxn * grd.nyn * grd.nzn);

    // implicit mover
    // #pragma omp parallel for // only if use mover_PC_V
    iMover = cpuSecond();  // start timer for mover

    for(int is=0; is<param.ns; is++){
        //cudaStreamSynchronize(streams[is]);
        setZeroDensities_gpu(&streams[is], &grd, grid_gpu_ptr, idn_gpu_ptr, ids_gpu_ptr[is], is);
    }

    for (int is=0; is < param.ns; is++)
    {
        //cudaStreamSynchronize(streams[is]);
        int b = batch_update_particles(
            &streams[is], 
            &part[is], 
            &part_positions_gpu[is], 
            part_positions_gpu_ptr[is], 
            part_info_gpu_ptr[is], 
            field_gpu_ptr,
            grid_gpu_ptr, 
            ids_gpu_ptr[is], 
            param_gpu_ptr, 
            batch_size
            );
        std::cout << "Move and interpolate on gpu. Species " << is << " - " << b << " batches " << std::endl;            

    }
    std::cout << "***********************" << std::endl;

    // Apply BCs
    for (int is = 0; is < param.ns; is++)
    {
        /**
         * synchronization needed, interpolation needs to finish before we can
         * apply boundary conditions
         */ 
        cudaStreamSynchronize(streams[is]);
        applyBCids_gpu(&streams[is], ids_gpu_ptr[is], grid_gpu_ptr, &grd, &param);
    }

    // sum over species
    cudaDeviceSynchronize();
    for (int is=0; is < param.ns; is++)
    {
        sumOverSpecies_gpu(idn_gpu_ptr, ids_gpu_ptr[is], grid_gpu_ptr, &grd);
        interp_DS_copy_to_host(&ids_gpu[is], &ids[is], &grd);
    }
    cudaDeviceSynchronize();
    interp_DNet_copy_to_host(&idn_gpu, &idn, &grd);
 
    dMover = cpuSecond() - iMover;
    eMover += dMover;  // stop timer for mover
    std::cout << "Mover and interpolation time: " << dMover << std::endl;

    // Maxwell solver
    iField = cpuSecond();  // start timer for the interpolation step

    // interpolate charge density from center to node
    applyBCscalarDensN(idn.rhon, &grd, &param);
    interpN2Cinterp(idn.rhoc, idn.rhon, &grd);
    // calculate hat functions rhoh and Jxh, Jyh, Jzh
    calculateHatDensities(&id_aux, &idn, ids, &field, &grd, &param);

    //  Poisson correction
    if (param.PoissonCorrection)
      divergenceCleaning(&grd, &field_aux, &field, &idn, &param);
    // Calculate E
    calculateE(&grd, &field_aux, &field, &id_aux, ids, &param);
    // Calculate B
    calculateB(&grd, &field_aux, &field, &param);
    dField = cpuSecond() - iField;
    eField += dField;  // stop timer for interpolation
    std::cout << "Solver time: " << dInterp << std::endl;

    // write E, B, rho to disk
    if (cycle % param.FieldOutputCycle == 0) {
      iIO = cpuSecond();
      VTK_Write_Vectors(cycle, &grd, &field);
      VTK_Write_Scalars(cycle, &grd, ids, &idn);
      eIO += (cpuSecond() - iIO);  // stop timer for interpolation
    }

  }  // end of one PIC cycle

  /// Release the resources
  // deallocate field
  grid_deallocate(&grd);
  field_deallocate(&grd, &field);
  // interp
  interp_dens_net_deallocate(&grd, &idn);
  interp_dens_aux_deallocate(&grd, &id_aux);

  // Deallocate interpolated densities and particles
  for (int is=0; is < param.ns; is++)
  {
      interp_dens_species_deallocate(&grd,&ids[is]);
      particle_deallocate_host(&part[is]);
  }

  std::cout << "DEALLOCATED CPU RESOURCES" << std::endl;

  // Deallocate GPU resources
  //streams
  for(int is=0; is<param.ns; is++){
      checkCudaErrors(cudaStreamDestroy(streams[is]));
  }

  // Parameters
  checkCudaErrors(cudaFree(param_gpu_ptr));
  // Grid
  grid_deallocate_device(&grid_gpu, grid_gpu_ptr);
  // EMfield
  field_deallocate_device(&field_gpu, field_gpu_ptr);

  // interpDensSpecies
  for(int is = 0; is < param.ns; is++)
  {
      interp_DS_dealloc_device(&ids_gpu[is], ids_gpu_ptr[is]);
      particles_info_dealloc_device(part_info_gpu_ptr[is]);
      particles_positions_dealloc_device(&part_positions_gpu[is], part_positions_gpu_ptr[is]);
  }
  interp_DNet_dealloc_device(&idn_gpu, idn_gpu_ptr);

  delete[] ids_gpu;
  delete[] ids_gpu_ptr;
  delete[] part_info_gpu;
  delete[] part_info_gpu_ptr;
  delete[] part_positions_gpu;
  delete[] part_positions_gpu_ptr;
  std::cout << "DEALLOCATED GPU RESOURCES" << std::endl;
 
  // stop timer
  double iElaps = cpuSecond() - iStart;

  // Print timing of simulation
  std::cout << std::endl;
  std::cout << "**************************************" << std::endl;
  std::cout << "   Tot. Simulation Time (s) = " << iElaps << std::endl;
  std::cout << "   Mover Time / Cycle   (s) = " << eMover / param.ncycles
            << std::endl;
  std::cout << "   Interp. Time / Cycle (s) = " << eInterp / param.ncycles
            << std::endl;
  std::cout << "   Field Time / Cycle   (s) = " << eField / param.ncycles
            << std::endl;
  std::cout << "   IO Time / Cycle      (s) = " << eIO / param.ncycles
            << std::endl;
  std::cout << "**************************************" << std::endl;

  // exit
  return 0;
}
