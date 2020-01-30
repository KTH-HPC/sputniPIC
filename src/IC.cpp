#include "IC.h"
#include <math.h>
#include <iostream>
#include <sys/stat.h>

/** initialize for magnetic reconnection probelm with Harris current sheet */
void initGEM(struct parameters *param, struct grid *grd, struct EMfield *field,
             struct EMfield_aux *field_aux, struct particles *part,
             struct interpDensSpecies *ids) {
  // perturbation localized in X
  double pertX = 0.4;
  double xpert, ypert, exp_pert;

  // print settings
  std::cout << "*************************************************" << std::endl;
  std::cout << "**  Initialize GEM Challenge with Pertubation  **" << std::endl;
  std::cout << "*************************************************" << std::endl;
  std::cout << "** B0x = " << param->B0x << std::endl;
  std::cout << "** B0y = " << param->B0y << std::endl;
  std::cout << "** B0z = " << param->B0z << std::endl;
  std::cout << "** Delta (current sheet thickness) = " << param->delta
            << std::endl;
  for (int is = 0; is < param->ns; is++) {
    std::cout << "** rho species " << is << " = " << param->rhoINIT[is];
    if (is < 2)
      std::cout << " CURRENT SHEET " << std::endl;
    else
      std::cout << " BACKGROUND " << std::endl;
  }
  std::cout << "*************************************************" << std::endl;
  struct stat sb;

  if (param->RESTART != true) {
    /////////////////////////////////////////////////
    //////   FIELD AND DENSITY
    /////////////////////////////////////////////////
    // Set the electric field, magnetic field + rhos
#pragma omp parallel for private(xpert,ypert,exp_pert)
    for (int i = 0; i < grd->nxn; i++)
      for (int j = 0; j < grd->nyn; j++)
        for (int k = 0; k < grd->nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < param->ns; is++) {
            if (is < 2)  // current sheet
              ids[is].rhon[i][j][k] =
                  (FPinterp)(
                      (param->rhoINIT[is] /
                       (cosh((grd->YN[i][j][k] - grd->Ly / 2) / param->delta) *
                        cosh((grd->YN[i][j][k] - grd->Ly / 2) / param->delta)))) /
                  param->fourpi;
            else  // background
              ids[is].rhon[i][j][k] =
                  (FPinterp)param->rhoINIT[is] / param->fourpi;
          }
          // std::cout << "OK" << std::endl;
          // electric field
          // electric field
          field->Ex[i][j][k] = 0.0;
          field->Ey[i][j][k] = 0.0;
          field->Ez[i][j][k] = 0.0;
          field_aux->Exth[i][j][k] = 0.0;
          field_aux->Eyth[i][j][k] = 0.0;
          field_aux->Ezth[i][j][k] = 0.0;
          // Magnetic field
          field->Bxn[i][j][k] =
              param->B0x * tanh((grd->YN[i][j][k] - grd->Ly / 2) / param->delta);
          // add the initial GEM perturbation
          // Bxn[i][j][k] +=
          // (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXN(i,j,k)/Lx)*sin(M_PI*(grid->getYN(i,j,k)-
          // Ly/2)/Ly  );
          field->Byn[i][j][k] =
              param
                  ->B0y;  // -
                          // (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXN(i,j,k)/Lx)*cos(M_PI*(grid->getYN(i,j,k)-
                          // Ly/2)/Ly);
          // add the initial X perturbation
          xpert = grd->XN[i][j][k] - grd->Lx / 2;
          ypert = grd->YN[i][j][k] - grd->Ly / 2;
          exp_pert = exp(-(xpert / param->delta) * (xpert / param->delta) -
                         (ypert / param->delta) * (ypert / param->delta));
          field->Bxn[i][j][k] +=
              (param->B0x * pertX) * exp_pert *
              (-cos(M_PI * xpert / 10.0 / param->delta) *
                   cos(M_PI * ypert / 10.0 / param->delta) * 2.0 * ypert /
                   param->delta -
               cos(M_PI * xpert / 10.0 / param->delta) *
                   sin(M_PI * ypert / 10.0 / param->delta) * M_PI / 10.0);
          field->Byn[i][j][k] +=
              (param->B0x * pertX) * exp_pert *
              (cos(M_PI * xpert / 10.0 / param->delta) *
                   cos(M_PI * ypert / 10.0 / param->delta) * 2.0 * xpert /
                   param->delta +
               sin(M_PI * xpert / 10.0 / param->delta) *
                   cos(M_PI * ypert / 10.0 / param->delta) * M_PI / 10.0);
          // guide field
          field->Bzn[i][j][k] = param->B0z;
        }
    // calculate B and rho at centers cells: first argument is on center cell
    interpN2Cfield(field_aux->Bxc, field_aux->Byc, field_aux->Bzc, field->Bxn,
                   field->Byn, field->Bzn, grd);
  
    // interpolate densities species from node
    for (int is = 0; is < param->ns; is++) {
      interpN2Crho(&ids[is], grd);
    }
    
    /////////////////////////////////////////////////
    //////   PARTICLE
    /////////////////////////////////////////////////
    
    double harvest;
    double prob, theta, sign;
    long long counter;
    unsigned int seed;
    
    // loop over the species
    for (int is = 0; is < param->ns; is++) {
      // set particle counter to zero
      counter = 0;
      // set the seed for random number generator equal to species id
      seed = is;
      srand(seed);
      for (int i = 1; i < grd->nxc - 1; i++)
        for (int j = 1; j < grd->nyc - 1; j++)
          for (int k = 1; k < grd->nzc - 1; k++)
            for (int ii = 0; ii < part[is].npcelx; ii++)
              for (int jj = 0; jj < part[is].npcely; jj++)
                for (int kk = 0; kk < part[is].npcely; kk++) {
                  // initialize each particle position and charge. Particle
                  // uniform in space
                  part[is].x[counter] =
                      (ii + .5) * (grd->dx / part[is].npcelx) + grd->XN[i][j][k];
                  part[is].y[counter] =
                      (jj + .5) * (grd->dy / part[is].npcely) + grd->YN[i][j][k];
                  part[is].z[counter] =
                      (kk + .5) * (grd->dz / part[is].npcelz) + grd->ZN[i][j][k];
                  // q = charge * statistical weight
                  part[is].q[counter] = (part[is].qom / fabs(part[is].qom)) *
                                        (ids[is].rhoc[i][j][k] / part[is].npcel) *
                                        (1.0 / grd->invVOL);
      
                  //////////////// Maxwellian ////////////////
                  // u
                  harvest = rand() / (double)RAND_MAX;
                  prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
                  harvest = rand() / (double)RAND_MAX;
                  theta = 2.0 * M_PI * harvest;
                  part[is].u[counter] =
                      part[is].u0 + part[is].uth * prob * cos(theta);
                  // check u
                  if (part[is].u[counter] > param->c) {
                    std::cout << "ERROR - u VELOCITY > c !" << std::endl;
                    exit(EXIT_FAILURE);
                  }
                  // u and w
                  harvest = rand() / (double)RAND_MAX;
                  prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
                  harvest = rand() / (double)RAND_MAX;
                  theta = 2.0 * M_PI * harvest;
                  part[is].w[counter] =
                      part[is].w0 + part[is].wth * prob * cos(theta);
                  part[is].v[counter] =
                      part[is].v0 + part[is].vth * prob * sin(theta);
                  // check v and w
                  if (part[is].v[counter] > param->c) {
                    std::cout << "ERROR - v VELOCITY > c !" << std::endl;
                    exit(EXIT_FAILURE);
                  }
                  if (part[is].w[counter] > param->c) {
                    std::cout << "ERROR - w VELOCITY > c !" << std::endl;
                    exit(EXIT_FAILURE);
                  }
                  //  update particle counter
                  counter++;
                }  // end of one particles initialization
      
      }  // end of species initialization
    if (stat(param->RestartDirName.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
      save_ic_data(field, field_aux, grd, ids, part, param);
    }
    else {
      std::cout << "folder " << param->RestartDirName << " does not exist, not recording IC." << std::endl;
    }
  }
  else {
    read_ic_data(field, field_aux, grd, ids, part, param);
  }
}

void save_ic_data(struct EMfield *field, struct EMfield_aux *field_aux, struct grid *grd, struct interpDensSpecies *ids, struct particles *part, struct parameters *param)
{
  std::cout << "Writing ic_data..." << std::endl;
  FILE * pFile;
  size_t ret;
  std::string file_name;
  /* load Bx By Bz and friends */
  file_name = param->RestartDirName + "/Bxn_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "wb");
  ret = fwrite(field->Bxn_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Byn_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "wb");
  ret = fwrite(field->Byn_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Bzn_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "wb");
  ret = fwrite(field->Bzn_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Bxc_" + std::to_string(grd->nxc) + "_" + std::to_string(grd->nyc) + "_" + std::to_string(grd->nzc) + ".bin";
  pFile = fopen(file_name.c_str(), "wb");
  ret = fwrite(field_aux->Bxc_flat, sizeof(FPfield), grd->nxc*grd->nyc*grd->nzc, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Byc_" + std::to_string(grd->nxc) + "_" + std::to_string(grd->nyc) + "_" + std::to_string(grd->nzc) + ".bin";
  pFile = fopen(file_name.c_str(), "wb");
  ret = fwrite(field_aux->Byc_flat, sizeof(FPfield), grd->nxc*grd->nyc*grd->nzc, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Bzc_" + std::to_string(grd->nxc) + "_" + std::to_string(grd->nyc) + "_" + std::to_string(grd->nzc) + ".bin";
  pFile = fopen(file_name.c_str(), "wb");
  ret = fwrite(field_aux->Bzc_flat, sizeof(FPfield), grd->nxc*grd->nyc*grd->nzc, pFile);
  fclose(pFile);

  /* load Ex Ey Ez and friends */
  file_name = param->RestartDirName + "/Exn_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "wb");
  ret = fwrite(field->Ex_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Eyn_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "wb");
  ret = fwrite(field->Ey_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Ezn_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "wb");
  ret = fwrite(field->Ez_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Exth_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "wb");
  ret = fwrite(field_aux->Exth_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Eyth_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "wb");
  ret = fwrite(field_aux->Eyth_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Ezth_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "wb");
  ret = fwrite(field_aux->Ezth_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  /* load particle (x,y,z) (u,v,w) q */
  for (int is = 0; is < param->ns; is++) {
    /* load ids rhon, rhoc */
    file_name = param->RestartDirName + "/rhon_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "wb");
    ret = fwrite(ids[is].rhon_flat, sizeof(FPinterp), grd->nxn*grd->nyn*grd->nzn, pFile);
    fclose(pFile);
  
    file_name = param->RestartDirName + "/rhoc_" + std::to_string(grd->nxc) + "_" + std::to_string(grd->nyc) + "_" + std::to_string(grd->nzc) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "wb");
    ret = fwrite(ids[is].rhoc_flat, sizeof(FPinterp), grd->nxc*grd->nyc*grd->nzc, pFile);
    fclose(pFile);

    /* x */
    file_name = param->RestartDirName + "/x_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "wb");
    ret = fwrite(part[is].x, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);

    /* y */
    file_name = param->RestartDirName + "/y_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "wb");
    ret = fwrite(part[is].y, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);

    /* z */
    file_name = param->RestartDirName + "/z_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "wb");
    ret = fwrite(part[is].z, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);

    /* u */
    file_name = param->RestartDirName + "/u_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "wb");
    ret = fwrite(part[is].u, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);

    /* v */
    file_name = param->RestartDirName + "/v_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "wb");
    ret = fwrite(part[is].v, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);

    /* w */
    file_name = param->RestartDirName + "/w_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "wb");
    ret = fwrite(part[is].w, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);

    /* q */
    file_name = param->RestartDirName + "/q_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "wb");
    ret = fwrite(part[is].q, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);
  }
}


void read_ic_data(struct EMfield *field, struct EMfield_aux *field_aux, struct grid *grd, struct interpDensSpecies *ids, struct particles *part, struct parameters *param)
{
  std::cout << "Reading ic_data..." << std::endl;
  FILE * pFile;
  size_t ret;
  std::string file_name;
  /* load Bx By Bz and friends */
  file_name = param->RestartDirName + "/Bxn_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "rb");
  ret = fread(field->Bxn_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Byn_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "rb");
  ret = fread(field->Byn_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Bzn_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "rb");
  ret = fread(field->Bzn_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Bxc_" + std::to_string(grd->nxc) + "_" + std::to_string(grd->nyc) + "_" + std::to_string(grd->nzc) + ".bin";
  pFile = fopen(file_name.c_str(), "rb");
  ret = fread(field_aux->Bxc_flat, sizeof(FPfield), grd->nxc*grd->nyc*grd->nzc, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Byc_" + std::to_string(grd->nxc) + "_" + std::to_string(grd->nyc) + "_" + std::to_string(grd->nzc) + ".bin";
  pFile = fopen(file_name.c_str(), "rb");
  ret = fread(field_aux->Byc_flat, sizeof(FPfield), grd->nxc*grd->nyc*grd->nzc, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Bzc_" + std::to_string(grd->nxc) + "_" + std::to_string(grd->nyc) + "_" + std::to_string(grd->nzc) + ".bin";
  pFile = fopen(file_name.c_str(), "rb");
  ret = fread(field_aux->Bzc_flat, sizeof(FPfield), grd->nxc*grd->nyc*grd->nzc, pFile);
  fclose(pFile);

  /* load Ex Ey Ez and friends */
  file_name = param->RestartDirName + "/Exn_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "rb");
  ret = fread(field->Ex_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Eyn_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "rb");
  ret = fread(field->Ey_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Ezn_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "rb");
  ret = fread(field->Ez_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Exth_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "rb");
  ret = fread(field_aux->Exth_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Eyth_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "rb");
  ret = fread(field_aux->Eyth_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  file_name = param->RestartDirName + "/Ezth_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + ".bin";
  pFile = fopen(file_name.c_str(), "rb");
  ret = fread(field_aux->Ezth_flat, sizeof(FPfield), grd->nxn*grd->nyn*grd->nzn, pFile);
  fclose(pFile);

  /* load particle (x,y,z) (u,v,w) q */
  for (int is = 0; is < param->ns; is++) {
    /* load ids rhon, rhoc */
    file_name = param->RestartDirName + "/rhon_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "rb");
    ret = fread(ids[is].rhon_flat, sizeof(FPinterp), grd->nxn*grd->nyn*grd->nzn, pFile);
    fclose(pFile);
  
    file_name = param->RestartDirName + "/rhoc_" + std::to_string(grd->nxc) + "_" + std::to_string(grd->nyc) + "_" + std::to_string(grd->nzc) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "rb");
    ret = fread(ids[is].rhoc_flat, sizeof(FPinterp), grd->nxc*grd->nyc*grd->nzc, pFile);
    fclose(pFile);

    /* x */
    file_name = param->RestartDirName + "/x_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "rb");
    ret = fread(part[is].x, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);

    /* y */
    file_name = param->RestartDirName + "/y_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "rb");
    ret = fread(part[is].y, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);

    /* z */
    file_name = param->RestartDirName + "/z_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "rb");
    ret = fread(part[is].z, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);

    /* u */
    file_name = param->RestartDirName + "/u_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "rb");
    ret = fread(part[is].u, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);

    /* v */
    file_name = param->RestartDirName + "/v_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "rb");
    ret = fread(part[is].v, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);

    /* w */
    file_name = param->RestartDirName + "/w_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "rb");
    ret = fread(part[is].w, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);

    /* q */
    file_name = param->RestartDirName + "/q_" + std::to_string(grd->nxn) + "_" + std::to_string(grd->nyn) + "_" + std::to_string(grd->nzn) + "_" + std::to_string(is) + "_" + std::to_string(part[is].npmax) + ".bin";
    pFile = fopen(file_name.c_str(), "rb");
    ret = fread(part[is].q, sizeof(FPpart), part[is].npmax, pFile);
    fclose(pFile);
  }
}

/** initialize uniform electric and magnetic field */
void initUniform(struct parameters *param, struct grid *grd,
             struct EMfield *field, struct EMfield_aux *field_aux,
             struct particles *part, struct interpDensSpecies *ids) {
// perturbation localized in X
double pertX = 0.4;
double xpert, ypert, exp_pert;

// print settings
std::cout << "*************************************************" << std::endl;
std::cout << "**  Initialize UNIFORM Plasma Distribution     **" << std::endl;
std::cout << "*************************************************" << std::endl;
std::cout << "** B0x = " << param->B0x << std::endl;
std::cout << "** B0y = " << param->B0y << std::endl;
std::cout << "** B0z = " << param->B0z << std::endl;
for (int is = 0; is < param->ns; is++) {
std::cout << "** rho species " << is << " = " << param->rhoINIT[is];
}
std::cout << "*************************************************" << std::endl;

/////////////////////////////////////////////////
//////   FIELD AND DENSITY
/////////////////////////////////////////////////
// Set the electric field, magnetic field + rhos
for (int i = 0; i < grd->nxn; i++)
for (int j = 0; j < grd->nyn; j++)
  for (int k = 0; k < grd->nzn; k++) {
    // initialize the density for species
    for (int is = 0; is < param->ns; is++) {
      ids[is].rhon[i][j][k] = (FPinterp)param->rhoINIT[is] / param->fourpi;
    }
    // std::cout << "OK" << std::endl;
    // electric field
    // electric field
    field->Ex[i][j][k] = 0.0;
        field->Ey[i][j][k] = 0.0;
        field->Ez[i][j][k] = 0.0;
        field_aux->Exth[i][j][k] = 0.0;
        field_aux->Eyth[i][j][k] = 0.0;
        field_aux->Ezth[i][j][k] = 0.0;
        // Magnetic field
        field->Bxn[i][j][k] = param->B0x;
        field->Byn[i][j][k] = param->B0y;
        field->Bzn[i][j][k] = param->B0z;
      }
  // calculate B and rho at centers cells: first argument is on center cell
  interpN2Cfield(field_aux->Bxc, field_aux->Byc, field_aux->Bzc, field->Bxn,
                 field->Byn, field->Bzn, grd);

  // interpolate densities species from node
  for (int is = 0; is < param->ns; is++) {
    interpN2Crho(&ids[is], grd);
  }

  /////////////////////////////////////////////////
  //////   PARTICLES
  /////////////////////////////////////////////////

  double harvest;
  double prob, theta, sign;
  long long counter;

  // loop over the species
  for (int is = 0; is < param->ns; is++) {
    // set particle counter to zero
    counter = 0;
    // set the seed for random number generator equal to species id
    srand(is);
    for (int i = 1; i < grd->nxc - 1; i++)
      for (int j = 1; j < grd->nyc - 1; j++)
        for (int k = 1; k < grd->nzc - 1; k++)
          for (int ii = 0; ii < part[is].npcelx; ii++)
            for (int jj = 0; jj < part[is].npcely; jj++)
              for (int kk = 0; kk < part[is].npcely; kk++) {
                // initialize each particle position and charge. Particle
                // uniform in space
                part[is].x[counter] =
                    (ii + .5) * (grd->dx / part[is].npcelx) + grd->XN[i][j][k];
                part[is].y[counter] =
                    (jj + .5) * (grd->dy / part[is].npcely) + grd->YN[i][j][k];
                part[is].z[counter] =
                    (kk + .5) * (grd->dz / part[is].npcelz) + grd->ZN[i][j][k];
                // q = charge * statistical weight
                part[is].q[counter] = (part[is].qom / fabs(part[is].qom)) *
                                      (ids[is].rhoc[i][j][k] / part[is].npcel) *
                                      (1.0 / grd->invVOL);

                //////////////// Maxwellian ////////////////
                // u
                harvest = rand() / (double)RAND_MAX;
                prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
                harvest = rand() / (double)RAND_MAX;
                theta = 2.0 * M_PI * harvest;
                part[is].u[counter] =
                    part[is].u0 + part[is].uth * prob * cos(theta);
                // check u
                if (part[is].u[counter] > param->c) {
                  std::cout << "ERROR - u VELOCITY > c !" << std::endl;
                  exit(EXIT_FAILURE);
                }
                // u and w
                harvest = rand() / (double)RAND_MAX;
                prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
                harvest = rand() / (double)RAND_MAX;
                theta = 2.0 * M_PI * harvest;
                part[is].w[counter] =
                    part[is].w0 + part[is].wth * prob * cos(theta);
                part[is].v[counter] =
                    part[is].v0 + part[is].vth * prob * sin(theta);
                // check v and w
                if (part[is].v[counter] > param->c) {
                  std::cout << "ERROR - v VELOCITY > c !" << std::endl;
                  exit(EXIT_FAILURE);
                }
                if (part[is].w[counter] > param->c) {
                  std::cout << "ERROR - w VELOCITY > c !" << std::endl;
                  exit(EXIT_FAILURE);
                }
                //  update particle counter
                counter++;
              }  // end of one particles initialization

  }  // end of species initialization
}
