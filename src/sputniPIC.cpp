/** A mixed-precision implicit Particle-in-Cell simulator for heterogeneous
 * systems **/

#include <omp.h>

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

#include "gpu/cuda_helper.h"

int main(int argc, char** argv) {
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
  double dMover = 0.0, dInterp = 0.0, dField = 0.0, dIO = 0.0;

  int num_devices;
  checkCudaErrors(cudaGetDeviceCount(&num_devices));

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
  particles* part = new particles[param.ns];
  // allocation
  for (int is = 0; is < param.ns; is++) {
    particle_allocate_host(&param, &part[is], is);
  }

  // Initialization
  initGEM(&param, &grd, &field, &field_aux, part, ids);
  // initUniform(&param,&grd,&field,&field_aux,part,ids);

  // ********************** GPU ALLOCATION ********************//
  // Create the parameters on the GPU and copy the data to the device.
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
    field_alloc_and_copy_to_device(&grd, &field_gpu[device_id], &field,
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

  int np = 0;
  for (int is = 0; is < param.ns; is++) {
    np += part[is].nop;
  }

  int batch_size = get_appropriate_batch_size(param.ns);
  if (batch_size <= 0) {
    return -1;
  }
  batch_size *= num_devices; /* size computed based on free mem of one device */

#ifdef MEMCHECK
  std::cout << "CUDA return check: Enabled" << std::endl;
#else
  std::cout << "CUDA return check: Enabled" << std::endl;
#endif
  std::cout << "Total number of cores: " << omp_get_max_threads() << std::endl;
  std::cout << "Total number of GPUs: " << num_devices << std::endl;
  std::cout << "Total number of particles: " << np << "; "
            << (np * sizeof(FPpart) * 6) / (1024 * 1024) << " MB of data."
            << std::endl;
  std::cout << "Allocating "
            << (batch_size * param.ns * sizeof(FPpart) * 6) / (1024 * 1024)
            << " MB of memory for particles on gpu" << std::endl;
  std::cout << "(batch_size of "
            << (batch_size * sizeof(FPpart) * 6) / (1024 * 1024) << " MB)"
            << std::endl;

  // allocate field and copy to GPUs
  for (int is = 0; is < param.ns; is++) {
    int device_id = (is + num_devices) % num_devices;
    checkCudaErrors(cudaSetDevice(device_id));
    particles_info_alloc_and_copy_to_device(&part_info_gpu[is],
                                            &part_info_gpu_ptr[is], &part[is]);
    particles_positions_alloc_device(&part_positions_gpu[is],
                                     &part_positions_gpu_ptr[is], batch_size);
  }

  // create cuda streams for asynchronous workloading
  cudaStream_t* streams = new cudaStream_t[param.ns];
  for (int is = 0; is < param.ns; is++) {
    int device_id = (is + num_devices) % num_devices;
    checkCudaErrors(cudaSetDevice(device_id));
    checkCudaErrors(
        cudaStreamCreateWithFlags(&streams[is], cudaStreamNonBlocking));
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
    setZeroNetDensities(&idn, &grd);

    // implicit mover
    iMover = cpuSecond();  // start timer for mover

    // blocking copy updated field from prev cycle
    for (int device_id = 0; device_id < num_devices; device_id++) {
      checkCudaErrors(cudaSetDevice(device_id));
      copy_from_host_to_device(&field_gpu[device_id], &field,
                               grd.nxn * grd.nyn * grd.nzn);
    }

    // async set zero, move and interp
    for (int is = 0; is < param.ns; is++) {
      int device_id = (is + num_devices) % num_devices;
      checkCudaErrors(cudaSetDevice(device_id));
      setZeroSpeciesDensities_gpu(&streams[is], &grd, grid_gpu_ptr[device_id],
                                  &ids_gpu[is], is);
      int b = batch_update_particles(
          &streams[is], &part[is], &part_positions_gpu[is],
          part_positions_gpu_ptr[is], part_info_gpu_ptr[is],
          field_gpu_ptr[device_id], grid_gpu_ptr[device_id], ids_gpu_ptr[is],
          param_gpu_ptr[device_id], batch_size);

      std::cout << "Move and interpolate on gpu " << device_id << " Species "
                << is << " - " << b << " batches " << std::endl;
    }
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

    dMover = cpuSecond() - iMover;
    eMover += dMover;  // stop timer for mover

    // sum over species
    sumOverSpecies(&idn, ids, &grd, param.ns);

    // interpolate charge density from center to node
    applyBCscalarDensN(idn.rhon, &grd, &param);
    interpN2Cinterp(idn.rhoc, idn.rhon, &grd);
    // calculate hat functions rhoh and Jxh, Jyh, Jzh
    calculateHatDensities(&id_aux, &idn, ids, &field, &grd, &param);

    // Start Maxwell solver
    iField = cpuSecond();  // start timer for the interpolation step

    //  Poisson correction
    if (param.PoissonCorrection)
      divergenceCleaning(&grd, &field_aux, &field, &idn, &param);
    // Calculate E
    calculateE(&grd, &field_aux, &field, &id_aux, ids, &param);
    // Calculate B
    calculateB(&grd, &field_aux, &field, &param);

    dField = cpuSecond() - iField;
    eField += dField;  // stop timer for solvers

    // write E, B, rho to disk
    if (cycle % param.FieldOutputCycle == 0) {
      iIO = cpuSecond();
      VTK_Write_Vectors(cycle, &grd, &field);
      VTK_Write_Scalars(cycle, &grd, ids, &idn);
      dIO = (cpuSecond() - iIO);
      eIO += dIO;  // stop timer for interpolation
    }

    // print dummy zero for interp
    std::cout << "Timing Cycle " << cycle << " : " << dMover << " 0.0 " << dField
              << " " << dIO << std::endl;

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
    particles_positions_dealloc_device(&part_positions_gpu[is],
                                       part_positions_gpu_ptr[is]);
  }

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
