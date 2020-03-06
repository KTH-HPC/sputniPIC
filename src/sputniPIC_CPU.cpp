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

void update_statistics(double *mean, double *variance, double new_value, long count)
{
  double delta = new_value - *mean;
  *mean += delta / (double)count;
  double delta2 = new_value - *mean;
  *variance += delta * delta2;
}


int main(int argc, char **argv){

    std::cout << "Total number of cores: " << omp_get_max_threads() << std::endl;
#ifdef USE_GPU
    std::cout << "Version: GPU" << std::endl;
#else
    std::cout << "Version: CPU" << std::endl;
#endif
    
    // Read the inputfile and fill the param structure
    parameters param;
    // Read the input file name from command line
    readInputFile(&param,argc,argv);
    printParameters(&param);
    saveParameters(&param);
    
    // Timing variables
    double iStart = cpuSecond();
    double iMover, iInterp, iField, iIO, eMover = 0.0, eInterp= 0.0, eField = 0.0, eIO=0.0;
    double dMover, dInterp, dField, dIO;
    double avg_mover = 0.0, avg_interp = 0.0, avg_field = 0.0, avg_IO = 0.0;
    double stddev_mover = 0.0, stddev_interp = 0.0, stddev_field = 0.0, stddev_IO = 0.0;
    
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
    // allocation
    for (int is=0; is < param.ns; is++)
        particle_allocate(&param,&part[is],is);
       
    
    // Initialization
    initGEM(&param,&grd,&field,&field_aux,part,ids);
    //initUniform(&param,&grd,&field,&field_aux,part,ids);
    
    // **********************************************************//
    // **** Start the Simulation!  Cycle index start from 1  *** //
    // **********************************************************//
    for (int cycle = param.first_cycle_n; cycle < (param.first_cycle_n + param.ncycles); cycle++) {
        dMover = 0.0; dInterp = 0.0; dField = 0.0; dIO = 0.0;
        
        std::cout << std::endl;
        std::cout << "***********************" << std::endl;
        std::cout << "   cycle = " << cycle << std::endl;
        std::cout << "***********************" << std::endl;
    
        // set to zero the densities - needed for interpolation
        setZeroDensities(&idn,ids,&grd,param.ns);
        
        
        
        // implicit mover
        iMover = cpuSecond(); // start timer for mover
        // #pragma omp parallel for // only if use mover_PC_V
        for (int is=0; is < param.ns; is++)
            mover_PC(&part[is],&field,&grd,&param);
            //mover_PC_V(&part[is],&field,&grd,&param);
            //mover_interp(&part[is], &field, &ids[is],&grd, &param);
        dMover = (cpuSecond() - iMover);
        eMover += dMover; // stop timer for mover
        update_statistics(&avg_mover, &stddev_mover, dMover, cycle);
        
        
        
        std::cout << "*** INTERPOLATION P2G ***" << std::endl;
        // interpolation particle to grid
        iInterp = cpuSecond(); // start timer for the interpolation step
        // interpolate species: MAXIMUM parallelism is number of species
        #pragma omp parallel for
        for (int is=0; is < param.ns; is++)
            interpP2G(&part[is],&ids[is],&grd);
        dInterp = (cpuSecond() - iInterp); // stop timer for interpolation
        eInterp += dInterp; // stop timer for interpolation
        update_statistics(&avg_interp, &stddev_interp, dInterp, cycle);

        // apply BC to interpolated densities
        for (int is=0; is < param.ns; is++)
            applyBCids(&ids[is],&grd,&param);
        // sum over species
        sumOverSpecies(&idn,ids,&grd,param.ns);
        // interpolate charge density from center to node
        applyBCscalarDensN(idn.rhon,&grd,&param);
        interpN2Cinterp(idn.rhoc,idn.rhon, &grd);
        // calculate hat functions rhoh and Jxh, Jyh, Jzh
        calculateHatDensities(&id_aux, &idn, ids, &field, &grd, &param);
        
        
        
        // Maxwell solver
        iField = cpuSecond(); // start timer for the interpolation step
        //  Poisson correction
        if (param.PoissonCorrection)
            divergenceCleaning(&grd,&field_aux,&field,&idn,&param);
        // Calculate E
        calculateE(&grd,&field_aux,&field,&id_aux,ids,&param);
        // Calculate B
        calculateB(&grd,&field_aux,&field,&param);
        dField = (cpuSecond() - iField); // stop timer for interpolation
        eField += dField; // stop timer for interpolation
        update_statistics(&avg_field, &stddev_field, dField, cycle);
        
        
        // write E, B, rho to disk
        if (cycle%param.FieldOutputCycle==0){
            iIO = cpuSecond();
            VTK_Write_Vectors(cycle, &grd,&field, &param);
            VTK_Write_Scalars(cycle, &grd,ids,&idn, &param);
            dIO = (cpuSecond() - iIO); // stop timer for interpolation
            eIO += dIO; // stop timer for interpolation
            update_statistics(&avg_IO, &stddev_IO, dIO, cycle);
        }

        std::cout << "Timing Cycle " << cycle << " : " << dMover << " " << dInterp << " " << dField << " " << dIO << std::endl;
    
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
    
    
    // stop timer
    double iElaps = cpuSecond() - iStart;
    
    // Print timing of simulation
    std::cout << std::endl;
    std::cout << "**************************************" << std::endl;
    std::cout << "   Tot. Simulation Time (s) = " << iElaps << std::endl;
    std::cout << "   Mover Time / Cycle   (s) = " << eMover/param.ncycles << std::endl;
    std::cout << "   Interp. Time / Cycle (s) = " << eInterp/param.ncycles  << std::endl;
    std::cout << "   Field Time / Cycle   (s) = " << eField/param.ncycles  << std::endl;
    std::cout << "   IO Time / Cycle      (s) = " << eIO/param.ncycles  << std::endl;
    std::cout << "**************************************" << std::endl;

    std::cout << "Mover: " << avg_mover << " " << sqrt(stddev_mover / (param.ncycles - 1)) << std::endl;
    std::cout << "Interp: " << avg_interp << " " << sqrt(stddev_interp / (param.ncycles - 1)) << std::endl;
    std::cout << "Field: " << avg_field << " " << sqrt(stddev_field / (param.ncycles - 1)) << std::endl;
    std::cout << "IO: " << avg_IO << " " << sqrt(stddev_IO / (param.ncycles -1)) << std::endl;

    
    // exit
    return 0;
}


