#include <string.h>
#include <cmath>

#include "ConfigFile.h"
#include "RW_IO.h"
#include "input_array.h"

/** read the inputfile given via the command line */
void readInputFile(struct parameters *param, int argc, char **argv) {
  // /////////////////////
  // Read the command line
  if (argc < 2) {
    std::cout << "Need to provide the input file for the sputniPIC simulation"
              << std::endl;
    exit(EXIT_FAILURE);
  } else if (argc < 3) {
    param->inputfile = argv[1];
    param->RESTART = false;
  } else {
    if (strcmp(argv[1], "restart") == 0) {
      param->inputfile = argv[2];
      param->RESTART = true;
    } else if (strcmp(argv[2], "restart") == 0) {
      param->inputfile = argv[1];
      param->RESTART = true;
    } else {
      std::cout << "Error when providing command line arguments" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // /////////////////////
  // Loading the input file
  ConfigFile config(param->inputfile);

  /** light speed */
  param->c = config.read<double>("c", 1);
  /** 4  pi */
  param->fourpi = config.read<double>("fourpi", 12.5663706144);
  /** time step */
  param->dt = config.read<double>("dt");
  /** decentering parameter */
  param->th = config.read<double>("th", 1);

  /** number of time cycles */
  param->ncycles = config.read<int>("ncycles");
  /** mover predictor correcto iteration */
  param->NiterMover = config.read<int>("NiterMover", 3);
  /** number of particle of subcycles in the mover */
  param->n_sub_cycles = config.read<int>("n_sub_cycles", 1);

  /** number of particle batches per species when run on GPU **/
  param->number_of_batches = config.read<int>("number_of_batches", 16);
  /** number of threads per block to use when running kernels on GPU **/
  param->threads_per_block = config.read<long>("threads_per_block", 256);

  /** simulation box length - X direction   */
  param->Lx = config.read<double>("Lx", 1);
  /** simulation box length - Y direction   */
  param->Ly = config.read<double>("Ly", 1);
  /** simulation box length - Z direction   */
  param->Lz = config.read<double>("Lz", 1);
  /** number of cells - X direction        */
  param->nxc = config.read<int>("nxc", 1);
  /** number of cells - Y direction        */
  param->nyc = config.read<int>("nyc", 1);
  /** number of cells - Z direction        */
  param->nzc = config.read<int>("nzc", 1);
  /** object center X, e.g. planet or comet   */
  param->x_center = config.read<double>("x_center", 0.5);
  /** object center Y, e.g. planet or comet   */
  param->y_center = config.read<double>("y_center", 0.5);
  /** object center Z, e.g. planet or comet   */
  param->z_center = config.read<double>("z_center", 0.5);
  /** object size - assuming a cubic box or sphere  */
  param->L_square = config.read<double>("L_square", 0.2);

  // THIS IS IMPORTANT
  /** number of actual species */
  param->ns = config.read<int>("ns");
  ////
  ////

  /** read input parameters for particles */

  // We read maximum 6 species from inputfile
  array_int npcelx0 = config.read<array_int>("npcelx");
  array_int npcely0 = config.read<array_int>("npcely");
  array_int npcelz0 = config.read<array_int>("npcelz");
  array_double qom0 = config.read<array_double>("qom");
  array_double uth0 = config.read<array_double>("uth");
  array_double vth0 = config.read<array_double>("vth");
  array_double wth0 = config.read<array_double>("wth");
  array_double u00 = config.read<array_double>("u0");
  array_double v00 = config.read<array_double>("v0");
  array_double w00 = config.read<array_double>("w0");

  param->npcelx[0] = npcelx0.a;
  param->npcely[0] = npcely0.a;
  param->npcelz[0] = npcelz0.a;
  param->qom[0] = qom0.a;
  param->uth[0] = uth0.a;
  param->vth[0] = vth0.a;
  param->wth[0] = wth0.a;
  param->u0[0] = u00.a;
  param->v0[0] = v00.a;
  param->w0[0] = w00.a;

  if (param->ns > 1) {
    param->npcelx[1] = npcelx0.b;
    param->npcely[1] = npcely0.b;
    param->npcelz[1] = npcelz0.b;
    param->qom[1] = qom0.b;
    param->uth[1] = uth0.b;
    param->vth[1] = vth0.b;
    param->wth[1] = wth0.b;
    param->u0[1] = u00.b;
    param->v0[1] = v00.b;
    param->w0[1] = w00.b;
  }
  if (param->ns > 2) {
    param->npcelx[2] = npcelx0.c;
    param->npcely[2] = npcely0.c;
    param->npcelz[2] = npcelz0.c;
    param->qom[2] = qom0.c;
    param->uth[2] = uth0.c;
    param->vth[2] = vth0.c;
    param->wth[2] = wth0.c;
    param->u0[2] = u00.c;
    param->v0[2] = v00.c;
    param->w0[2] = w00.c;
  }
  if (param->ns > 3) {
    param->npcelx[3] = npcelx0.d;
    param->npcely[3] = npcely0.d;
    param->npcelz[3] = npcelz0.d;
    param->qom[3] = qom0.d;
    param->uth[3] = uth0.d;
    param->vth[3] = vth0.d;
    param->wth[3] = wth0.d;
    param->u0[3] = u00.d;
    param->v0[3] = v00.d;
    param->w0[3] = w00.d;
  }
  if (param->ns > 4) {
    param->npcelx[4] = npcelx0.e;
    param->npcely[4] = npcely0.e;
    param->npcelz[4] = npcelz0.e;
    param->qom[4] = qom0.e;
    param->uth[4] = uth0.e;
    param->vth[4] = vth0.e;
    param->wth[4] = wth0.e;
    param->u0[4] = u00.e;
    param->v0[4] = v00.e;
    param->w0[4] = w00.e;
  }
  if (param->ns > 5) {
    param->npcelx[5] = npcelx0.f;
    param->npcely[5] = npcely0.f;
    param->npcelz[5] = npcelz0.f;
    param->qom[5] = qom0.f;
    param->uth[5] = uth0.f;
    param->vth[5] = vth0.f;
    param->wth[5] = wth0.f;
    param->u0[5] = u00.f;
    param->v0[5] = v00.f;
    param->w0[5] = w00.f;
  }

  // Initialization of densities
  array_double rhoINIT0 = config.read<array_double>("rhoINIT");
  param->rhoINIT[0] = rhoINIT0.a;
  if (param->ns > 1) param->rhoINIT[1] = rhoINIT0.b;
  if (param->ns > 2) param->rhoINIT[2] = rhoINIT0.c;
  if (param->ns > 3) param->rhoINIT[3] = rhoINIT0.d;
  if (param->ns > 4) param->rhoINIT[4] = rhoINIT0.e;
  if (param->ns > 5) param->rhoINIT[5] = rhoINIT0.f;

  // Calculate the total number of particles in the domain
  param->NpMaxNpRatio = config.read<double>("NpMaxNpRatio", 1.0);
  int npcel = 0;
  for (int i = 0; i < param->ns; i++) {
    npcel = param->npcelx[i] * param->npcely[i] * param->npcelz[i];
    param->np[i] = npcel * param->nxc * param->nyc * param->nzc;
    param->npMax[i] = (long)(param->NpMaxNpRatio * param->np[i]);
  }

  // Boundary Conditions
  /** Periodicity for fields X **/
  param->PERIODICX = config.read<bool>("PERIODICX", true);
  /** Periodicity for fields Y **/
  param->PERIODICY = config.read<bool>("PERIODICY", false);
  /** Periodicity for fields Z **/
  param->PERIODICZ = config.read<bool>("PERIODICZ", true);
  /** Periodicity for Particles X **/
  param->PERIODICX_P = config.read<bool>("PERIODICX_P", true);
  /** Periodicity for Particles Y **/
  param->PERIODICY_P = config.read<bool>("PERIODICY_P", false);
  /** Periodicity for Particles Y **/
  param->PERIODICZ_P = config.read<bool>("PERIODICZ_P", true);

  // PHI Electrostatic Potential
  param->bcPHIfaceXright = config.read<int>("bcPHIfaceXright", 1);
  param->bcPHIfaceXleft = config.read<int>("bcPHIfaceXleft", 1);
  param->bcPHIfaceYright = config.read<int>("bcPHIfaceYright", 1);
  param->bcPHIfaceYleft = config.read<int>("bcPHIfaceYleft", 1);
  param->bcPHIfaceZright = config.read<int>("bcPHIfaceZright", 1);
  param->bcPHIfaceZleft = config.read<int>("bcPHIfaceZleft", 1);

  // EM field boundary condition
  param->bcEMfaceXright = config.read<int>("bcEMfaceXright", 1);
  param->bcEMfaceXleft = config.read<int>("bcEMfaceXleft", 1);
  param->bcEMfaceYright = config.read<int>("bcEMfaceYright", 1);
  param->bcEMfaceYleft = config.read<int>("bcEMfaceYleft", 1);
  param->bcEMfaceZright = config.read<int>("bcEMfaceZright", 1);
  param->bcEMfaceZleft = config.read<int>("bcEMfaceZleft", 1);

  // Particles Boundary condition
  param->bcPfaceXright = config.read<int>("bcPfaceXright", 1);
  param->bcPfaceXleft = config.read<int>("bcPfaceXleft", 1);
  param->bcPfaceYright = config.read<int>("bcPfaceYright", 1);
  param->bcPfaceYleft = config.read<int>("bcPfaceYleft", 1);
  param->bcPfaceZright = config.read<int>("bcPfaceZright", 1);
  param->bcPfaceZleft = config.read<int>("bcPfaceZleft", 1);

  // take the injection of the particless
  param->Vinj = config.read<double>("Vinj", 0.0);

  // Initialization
  param->B0x = config.read<double>("B0x", 0.0);
  param->B0y = config.read<double>("B0y", 0.0);
  param->B0z = config.read<double>("B0z", 0.0);
  param->delta = config.read<double>("delta", 0.5);

  /** Smoothing quantities */
  param->SmoothON = config.read<bool>("SmoothON",
                                      true);  // Smoothing is ON by default
  /** Smoothing value*/
  param->SmoothValue = config.read<double>(
      "SmoothValue", 0.5);  // between 0 and 1, typically 0.5
  /** Ntimes: smoothing is applied */
  param->SmoothTimes = config.read<int>("SmoothTimes", 6);

  // Waves info
  param->Nwaves = config.read<int>("Nwaves", 1);
  param->dBoB0 = config.read<double>("dBoB0", 0.0);
  param->WaveFile = config.read<string>("WaveFile", "WaveFile.txt");
  param->energy = config.read<double>("energy", 0.018199864696222);
  param->pitch_angle = config.read<double>(
      "pitch_angle", 0.698131700797732);  // 40 degree default

  param->verbose = config.read<bool>("verbose", true);

  // Poisson Correction
  param->PoissonCorrection = config.read<bool>("PoissonCorrection", true);
  param->CGtol = config.read<double>("CGtol", 1E-3);
  param->GMREStol = config.read<double>("GMREStol", 1E-3);

  // needed for restart (in this case no restart)
  param->first_cycle_n = 1;

  // take the output cycles
  param->FieldOutputCycle = config.read<int>("FieldOutputCycle", 10);
  param->ParticlesOutputCycle = config.read<int>("ParticlesOutputCycle", 100);
  param->RestartOutputCycle = config.read<int>("RestartOutputCycle", 100000);
  param->DiagnosticsOutputCycle =
      config.read<int>("DiagnosticsOutputCycle", 10);

  param->SaveDirName = config.read<string>("SaveDirName");
  param->RestartDirName = config.read<string>("RestartDirName");
}

/** Print Simulation Parameters */
void printParameters(struct parameters *param) {
  std::cout << std::endl;
  std::cout << "-------------------------" << std::endl;
  std::cout << "sputniPIC Sim. Parameters" << std::endl;
  std::cout << "-------------------------" << std::endl;
  std::cout << "Number of species    = " << param->ns << std::endl;
  for (int i = 0; i < param->ns; i++)
    std::cout << "Number of particles of species " << i << " = " << param->np[i]
              << "\t (MAX = " << param->npMax[i] << ")"
              << "  QOM = " << param->qom[i] << std::endl;
  std::cout << "x-Length                 = " << param->Lx << std::endl;
  std::cout << "y-Length                 = " << param->Ly << std::endl;
  std::cout << "z-Length                 = " << param->Lz << std::endl;
  std::cout << "Number of cells (x)      = " << param->nxc << std::endl;
  std::cout << "Number of cells (y)      = " << param->nyc << std::endl;
  std::cout << "Number of cells (z)      = " << param->nzc << std::endl;
  std::cout << "Time step                = " << param->dt << std::endl;
  std::cout << "Number of cycles         = " << param->ncycles << std::endl;
  std::cout << "Results saved in: " << param->SaveDirName << std::endl;
}

/** Save Simulation Parameters */
void saveParameters(struct parameters *param) {
  string temp;
  temp = param->SaveDirName + "/sputniPICparameters.txt";

  std::ofstream my_file(temp.c_str());

  my_file << "-----------------------------" << std::endl;
  my_file << "- sputniPIC Sim. Parameters -" << std::endl;
  my_file << "-----------------------------" << std::endl;

  my_file << std::endl;

  my_file << "Number of species    = " << param->ns << std::endl;
  for (int i = 0; i < param->ns; i++)
    my_file << "Number of particles of species " << i << " = " << param->np[i]
            << "\t (MAX = " << param->npMax[i] << ")"
            << "  QOM = " << param->qom[i] << std::endl;
  my_file << "---------------------------" << std::endl;
  my_file << "x-Length                 = " << param->Lx << std::endl;
  my_file << "y-Length                 = " << param->Ly << std::endl;
  my_file << "z-Length                 = " << param->Lz << std::endl;
  my_file << "Number of cells (x)      = " << param->nxc << std::endl;
  my_file << "Number of cells (y)      = " << param->nyc << std::endl;
  my_file << "Number of cells (z)      = " << param->nzc << std::endl;
  my_file << "---------------------------" << std::endl;
  my_file << "Time step                = " << param->dt << std::endl;
  my_file << "Number of cycles         = " << param->ncycles << std::endl;
  my_file << "---------------------------" << std::endl;
  for (int is = 0; is < param->ns; is++)
    my_file << "rho init species " << is << " = " << param->rhoINIT[is]
            << std::endl;
  my_file << "current sheet thickness  = " << param->delta << std::endl;
  my_file << "B0x                      = " << param->B0x << std::endl;
  my_file << "BOy                      = " << param->B0y << std::endl;
  my_file << "B0z                      = " << param->B0z << std::endl;
  my_file << "---------------------------" << std::endl;
  my_file << "Smooth ON                = " << param->SmoothON << std::endl;
  my_file << "GMRES error tolerance    = " << param->GMREStol << std::endl;
  my_file << "CG error tolerance       = " << param->CGtol << std::endl;
  my_file << "Mover error tolerance    = " << param->NiterMover << std::endl;
  my_file << "---------------------------" << std::endl;
  my_file << "Results saved in: " << param->SaveDirName << std::endl;
  my_file << "Restart saved in: " << param->RestartDirName << std::endl;
  my_file << "---------------------" << std::endl;

  my_file.close();
}

void VTK_Write_Vectors(int cycle, struct grid *grd, struct EMfield *field, struct parameters *param) {
  // stream file to be opened and managed
  string filename = "E";
  string temp;
  std::stringstream cc;
  cc << cycle;
  temp = param->SaveDirName + "/" + filename + "_" + cc.str();
  temp += ".vtk";
  std::cout << "Opening file: " << temp << std::endl;

  int nxn = grd->nxn;
  int nyn = grd->nyn;
  int nzn = grd->nzn;

  double dx = grd->dx;
  double dy = grd->dy;
  double dz = grd->dz;

  std::ofstream my_fileE(temp.c_str());
  my_fileE << "# vtk DataFile Version 1.0" << std::endl;
  my_fileE << "E field" << std::endl;
  my_fileE << "ASCII" << std::endl;
  my_fileE << "DATASET STRUCTURED_POINTS" << std::endl;
  my_fileE << "DIMENSIONS " << (nxn - 3) << " " << (nyn - 3) << " " << (nzn - 3)
           << std::endl;
  my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
  my_fileE << "SPACING " << dx << " " << dy << " " << dz << std::endl;
  my_fileE << "POINT_DATA " << (nxn - 3) * (nyn - 3) * (nzn - 3) << std::endl;
  my_fileE << "VECTORS E float" << std::endl;

  double Ex = 0, Ey = 0, Ez = 0;

  for (int k = 1; k < nzn - 2; k++)
    for (int j = 1; j < nyn - 2; j++)
      for (int i = 1; i < nxn - 2; i++) {
        Ex = field->Ex[i][j][k];
        if (fabs(Ex) < 1E-8) Ex = 0.0;
        Ey = field->Ey[i][j][k];
        if (fabs(Ey) < 1E-8) Ey = 0.0;
        Ez = field->Ez[i][j][k];
        if (fabs(Ez) < 1E-8) Ez = 0.0;
        my_fileE << Ex << " " << Ey << " " << Ez << std::endl;
      }

  my_fileE.close();

  filename = "B";
  temp = param->SaveDirName + "/" + filename + "_" + cc.str();
  temp += ".vtk";
  std::cout << "Opening file: " << temp << std::endl;
  std::ofstream my_file2(temp.c_str());
  my_file2 << "# vtk DataFile Version 1.0" << std::endl;
  my_file2 << "B field" << std::endl;
  my_file2 << "ASCII" << std::endl;
  my_file2 << "DATASET STRUCTURED_POINTS" << std::endl;
  my_file2 << "DIMENSIONS " << (nxn - 3) << " " << (nyn - 3) << " " << (nzn - 3)
           << std::endl;
  my_file2 << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
  my_file2 << "SPACING " << dx << " " << dy << " " << dz << std::endl;
  my_file2 << "POINT_DATA " << (nxn - 3) * (nyn - 3) * (nzn - 3) << std::endl;
  my_file2 << "VECTORS B float" << std::endl;

  double Bx = 0, By = 0, Bz = 0;

  for (int k = 1; k < nzn - 2; k++)
    for (int j = 1; j < nyn - 2; j++)
      for (int i = 1; i < nxn - 2; i++) {
        Bx = field->Bxn[i][j][k];
        if (fabs(Bx) < 1E-8) Bx = 0.0;
        By = field->Byn[i][j][k];
        if (fabs(By) < 1E-8) By = 0.0;
        Bz = field->Bzn[i][j][k];
        if (fabs(Bz) < 1E-8) Bz = 0.0;
        my_file2 << Bx << " " << By << " " << Bz << std::endl;
      }

  my_file2.close();
}

void VTK_Write_Scalars(int cycle, struct grid *grd,
                       struct interpDensSpecies *ids,
                       struct interpDensNet *idn,
                       struct parameters *param) {
  // stream file to be opened and managed
  string filename = "rhoe";
  string temp;
  std::stringstream cc;
  cc << cycle;
  temp = param->SaveDirName + "/" + filename + "_" + cc.str();
  temp += ".vtk";
  std::cout << "Opening file: " << temp << std::endl;

  // get the number of nodes
  int nxn = grd->nxn;
  int nyn = grd->nyn;
  int nzn = grd->nzn;

  // get the grid spacing
  double dx = grd->dx;
  double dy = grd->dy;
  double dz = grd->dz;

  std::ofstream my_file(temp.c_str());
  my_file << "# vtk DataFile Version 1.0" << std::endl;
  my_file << "Electron Density - Current Sheet " << std::endl;
  my_file << "ASCII" << std::endl;
  my_file << "DATASET STRUCTURED_POINTS" << std::endl;
  my_file << "DIMENSIONS " << (nxn - 3) << " " << (nyn - 3) << " " << (nzn - 3)
          << std::endl;
  my_file << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
  my_file << "SPACING " << dx << " " << dy << " " << dz << std::endl;
  my_file << "POINT_DATA " << (nxn - 3) * (nyn - 3) * (nzn - 3) << std::endl;
  my_file << "SCALARS rhoe float" << std::endl;
  my_file << "LOOKUP_TABLE default" << std::endl;

  for (int k = 1; k < nzn - 2; k++)
    for (int j = 1; j < nyn - 2; j++)
      for (int i = 1; i < nxn - 2; i++) {
        my_file << ids[0].rhon[i][j][k] << std::endl;
      }

  my_file.close();

  filename = "rhoi";
  temp = param->SaveDirName +"/" + filename + "_" + cc.str();
  temp += ".vtk";
  std::cout << "Opening file: " << temp << std::endl;
  std::ofstream my_file2(temp.c_str());
  my_file2 << "# vtk DataFile Version 1.0" << std::endl;
  my_file2 << "Ion Density - Current Sheet" << std::endl;
  my_file2 << "ASCII" << std::endl;
  my_file2 << "DATASET STRUCTURED_POINTS" << std::endl;
  my_file2 << "DIMENSIONS " << (nxn - 3) << " " << (nyn - 3) << " " << (nzn - 3)
           << std::endl;
  my_file2 << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
  my_file2 << "SPACING " << dx << " " << dy << " " << dz << std::endl;
  my_file2 << "POINT_DATA " << (nxn - 3) * (nyn - 3) * (nzn - 3) << std::endl;
  my_file2 << "SCALARS rhoi float" << std::endl;
  my_file2 << "LOOKUP_TABLE default" << std::endl;

  for (int k = 1; k < nzn - 2; k++)
    for (int j = 1; j < nyn - 2; j++)
      for (int i = 1; i < nxn - 2; i++) {
        my_file2 << ids[1].rhon[i][j][k] << std::endl;
      }

  my_file2.close();

  filename = "rho_net";
  temp = param->SaveDirName + "/" + filename + "_" + cc.str();
  temp += ".vtk";
  std::cout << "Opening file: " << temp << std::endl;
  std::ofstream my_file1(temp.c_str());
  my_file1 << "# vtk DataFile Version 1.0" << std::endl;
  my_file1 << "Net Charge Density" << std::endl;
  my_file1 << "ASCII" << std::endl;
  my_file1 << "DATASET STRUCTURED_POINTS" << std::endl;
  my_file1 << "DIMENSIONS " << (nxn - 3) << " " << (nyn - 3) << " " << (nzn - 3)
           << std::endl;
  my_file1 << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
  my_file1 << "SPACING " << dx << " " << dy << " " << dz << std::endl;
  my_file1 << "POINT_DATA " << (nxn - 3) * (nyn - 3) * (nzn - 3) << std::endl;
  my_file1 << "SCALARS rhonet float" << std::endl;
  my_file1 << "LOOKUP_TABLE default" << std::endl;

  for (int k = 1; k < nzn - 2; k++)
    for (int j = 1; j < nyn - 2; j++)
      for (int i = 1; i < nxn - 2; i++) {
        my_file1 << idn->rhon[i][j][k] << std::endl;
      }

  my_file1.close();
}
