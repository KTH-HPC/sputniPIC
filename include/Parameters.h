#ifndef PARAMETERS_H
#define PARAMETERS_H

// standard libraries
#include <string>


#define NS_MAX 6

/** Simulation Parameters */
struct parameters {
    
    /** light speed */
    double c;
    /** 4  pi */
    double fourpi;
    /** time step */
    double dt;
    /** decentering parameter */
    double th;
    
    /** number of time cycles */
    int ncycles;
    /** mover predictor correcto iteration */
    int NiterMover;
    /** number of particle of subcycles in the mover */
    int n_sub_cycles;
    
    /** simulation box length - X direction   */
    double Lx;
    /** simulation box length - Y direction   */
    double Ly;
    /** simulation box length - Z direction   */
    double Lz;
    /** number of cells - X direction        */
    int nxc;
    /** number of cells - Y direction        */
    int nyc;
    /** number of cells - Z direction        */
    int nzc;
    /** object center X, e.g. planet or comet   */
    double x_center;
    /** object center Y, e.g. planet or comet   */
    double y_center;
    /** object center Z, e.g. planet or comet   */
    double z_center;
    /** object size - assuming a cubic box or sphere  */
    double L_square;
    
    
    
    /** number of actual species */
    int ns;
    
    // This for maximum NS_MAX species. To have more increase the array size in NS_MAX
    /** number of particles per cell - X direction */
    int npcelx[NS_MAX];
    /** number of particles per cell - Y direction */
    int npcely[NS_MAX];
    /** number of particles per cell - Z direction */
    int npcelz[NS_MAX];
    /** number of particles array for different species */
    long np[NS_MAX];
    /** maximum number of particles array for different species */
    long npMax[NS_MAX];
    /** max number of particles */
    double NpMaxNpRatio;
    /** charge to mass ratio array for different species */
    double qom[NS_MAX];
    /** charge to mass ratio array for different species */
    double rhoINIT[NS_MAX];
    /** thermal velocity  - Direction X  */
    double uth[NS_MAX];
    /** thermal velocity  - Direction Y  */
    double vth[NS_MAX];
    /** thermal velocity  - Direction Z  */
    double wth[NS_MAX];
    /** Drift velocity - Direction X     */
    double u0[NS_MAX];
    /** Drift velocity - Direction Y    */
    double v0[NS_MAX];
    /** Drift velocity - Direction Z     */
    double w0[NS_MAX];
    
    
    
    /** Boundary Condition: Periodicity **/
    // here you have to set the topology for the fields
    /** Periodicity for fields X **/
    bool PERIODICX;
    /** Periodicity for fields Y **/
    bool PERIODICY;
    /** Periodicity for fields Z **/
    bool PERIODICZ;
    /** Periodicity for Particles X **/
    bool PERIODICX_P;
    /** Periodicity for Particles Y **/
    bool PERIODICY_P;
    /** Periodicity for Particles Y **/
    bool PERIODICZ_P;
    
    
    /** Boundary condition on particles
     0 = exit
     1 = perfect mirror
     2 = riemission
     */
    /** Boundary Condition Particles: FaceXright */
    int bcPfaceXright;
    /** Boundary Condition Particles: FaceXleft */
    int bcPfaceXleft;
    /** Boundary Condition Particles: FaceYright */
    int bcPfaceYright;
    /** Boundary Condition Particles: FaceYleft */
    int bcPfaceYleft;
    /** Boundary Condition Particles: FaceYright */
    int bcPfaceZright;
    /** Boundary Condition Particles: FaceYleft */
    int bcPfaceZleft;
    
    
    /** Field Boundary Condition
     0 =  Dirichlet Boundary Condition: specifies the valueto take pn the boundary of the domain
     1 =  Neumann Boundary Condition: specifies the value of derivative to take on the boundary of the domain
     2 =  Periodic Condition
     */
    /** Boundary Condition Electrostatic Potential: FaceXright */
    int bcPHIfaceXright;
    /** Boundary Condition Electrostatic Potential:FaceXleft */
    int bcPHIfaceXleft;
    /** Boundary Condition Electrostatic Potential:FaceYright */
    int bcPHIfaceYright;
    /** Boundary Condition Electrostatic Potential:FaceYleft */
    int bcPHIfaceYleft;
    /** Boundary Condition Electrostatic Potential:FaceZright */
    int bcPHIfaceZright;
    /** Boundary Condition Electrostatic Potential:FaceZleft */
    int bcPHIfaceZleft;
    
    /** Boundary Condition EM Field: FaceXright */
    int bcEMfaceXright;
    /** Boundary Condition EM Field: FaceXleft */
    int bcEMfaceXleft;
    /** Boundary Condition EM Field: FaceYright */
    int bcEMfaceYright;
    /** Boundary Condition EM Field: FaceYleft */
    int bcEMfaceYleft;
    /** Boundary Condition EM Field: FaceZright */
    int bcEMfaceZright;
    /** Boundary Condition EM Field: FaceZleft */
    int bcEMfaceZleft;
    
    /** velocity of the injection from the wall */
    double Vinj;
    
    
    /** Initial Condition*/
    /** current sheet thickness */
    double delta;
    /* Initial B0x */
    double B0x;
    /* Initial B0y */
    double B0y;
    /* Initial B0y */
    double B0z;
    /** Number of waves present in the system: used for turbulence studies */
    int Nwaves;
    /** Perturbation amplitude: used for turbulence studies */
    double dBoB0;
    /** pitch angle */
    double pitch_angle;
    /** energy of the particle */
    double energy;
    
    /** Smoothing quantities */
    bool SmoothON;
    /** Smoothing value*/
    double SmoothValue; // between 0 and 1, typically 0.5
    /** Ntimes: smoothing is applied */
    int SmoothTimes;
    
    
    /** boolean value for verbose results */
    bool verbose;
    /** RESTART */
    bool RESTART;
    
    
    /** Poisson Correction */
    bool PoissonCorrection;
    /** CG solver stopping criterium tolerance */
    double CGtol;
    /** GMRES solver stopping criterium tolerance */
    double GMREStol;
    
    /** needed if restart */
    int first_cycle_n;
    
    /** Output for field */
    int FieldOutputCycle;
    /** Output for particles */
    int ParticlesOutputCycle;
    /** restart cycle */
    int RestartOutputCycle;
    /** Output for diagnostics */
    int DiagnosticsOutputCycle;
    
    /** inputfile */
    std::string inputfile;
    /** SaveDirName     */
    std::string SaveDirName;
    /** RestartDirName     */
    std::string RestartDirName;
    /** name of the file with wave amplitude and phases */
    std::string WaveFile;
    
};
#endif
