#include "MaxwellSolver.h"
#include "Alloc.h"
#include "Basic.h"

void MaxwellSource(FPfield *bkrylov, grid *grd, EMfield *field,
                   EMfield_aux *field_aux, interpDens_aux *id_aux,
                   parameters *param) {
  // get drid point nodes
  int nxn = grd->nxn;
  int nyn = grd->nyn;
  int nzn = grd->nzn;

  FPfield delt = param->c * param->th * param->dt;

  // temporary arrays
  FPfield ***tempX = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***tempY = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***tempZ = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***tempXN = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***tempYN = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***tempZN = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***temp2X = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***temp2Y = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***temp2Z = newArr3(FPfield, nxn, nyn, nzn);

  // set to zero
  eqValue(0.0, tempX, nxn, nyn, nzn);
  eqValue(0.0, tempY, nxn, nyn, nzn);
  eqValue(0.0, tempZ, nxn, nyn, nzn);
  eqValue(0.0, tempXN, nxn, nyn, nzn);
  eqValue(0.0, tempYN, nxn, nyn, nzn);
  eqValue(0.0, tempZN, nxn, nyn, nzn);
  eqValue(0.0, temp2X, nxn, nyn, nzn);
  eqValue(0.0, temp2Y, nxn, nyn, nzn);
  eqValue(0.0, temp2Z, nxn, nyn, nzn);
  // communicate

  // fixBforcefree(grid,vct);
  // fixBgem(grid, vct);

  // prepare curl of B for known term of Maxwell solver: for the source term
  applyBCscalarFieldC(field_aux->Bxc, grd, param);  // tangential Neumann
  applyBCscalarFieldC(field_aux->Byc, grd, param);  // normal is zero
  applyBCscalarFieldC(field_aux->Bzc, grd, param);  // tangential Neumann

  curlC2N(tempXN, tempYN, tempZN, field_aux->Bxc, field_aux->Byc,
          field_aux->Bzc, grd);

  scale(temp2X, id_aux->Jxh, -param->fourpi / param->c, nxn, nyn, nzn);
  scale(temp2Y, id_aux->Jyh, -param->fourpi / param->c, nxn, nyn, nzn);
  scale(temp2Z, id_aux->Jzh, -param->fourpi / param->c, nxn, nyn, nzn);

  sum(temp2X, tempXN, nxn, nyn, nzn);
  sum(temp2Y, tempYN, nxn, nyn, nzn);
  sum(temp2Z, tempZN, nxn, nyn, nzn);
  scale(temp2X, delt, nxn, nyn, nzn);
  scale(temp2Y, delt, nxn, nyn, nzn);
  scale(temp2Z, delt, nxn, nyn, nzn);

  // communicateCenterBC_P(nxc, nyc, nzc, rhoh, 2, 2, 2, 2, 2, 2, vct);
  applyBCscalarDensC(id_aux->rhoh, grd, param);
  gradC2N(tempX, tempY, tempZ, id_aux->rhoh, grd);

  scale(tempX, -delt * delt * param->fourpi, nxn, nyn, nzn);
  scale(tempY, -delt * delt * param->fourpi, nxn, nyn, nzn);
  scale(tempZ, -delt * delt * param->fourpi, nxn, nyn, nzn);

  // sum E, past values
  sum(tempX, field->Ex, nxn, nyn, nzn);
  sum(tempY, field->Ey, nxn, nyn, nzn);
  sum(tempZ, field->Ez, nxn, nyn, nzn);

  // sum curl(B) + jhat part
  sum(tempX, temp2X, nxn, nyn, nzn);
  sum(tempY, temp2Y, nxn, nyn, nzn);
  sum(tempZ, temp2Z, nxn, nyn, nzn);

  // physical space -> Krylov space
  phys2solver(bkrylov, tempX, tempY, tempZ, nxn, nyn, nzn);

  delArr3(tempX, nxn, nyn);
  delArr3(tempY, nxn, nyn);
  delArr3(tempZ, nxn, nyn);

  delArr3(tempXN, nxn, nyn);
  delArr3(tempYN, nxn, nyn);
  delArr3(tempZN, nxn, nyn);

  delArr3(temp2X, nxn, nyn);
  delArr3(temp2Y, nxn, nyn);
  delArr3(temp2Z, nxn, nyn);
}

/** Maxwell Image */
void MaxwellImage(FPfield *im, FPfield *vector, EMfield *field,
                  interpDensSpecies *ids, grid *grd, parameters *param) {
  FPfield beta, edotb, omcx, omcy, omcz, denom;
  FPfield delt = param->c * param->th * param->dt;

  // get drid point nodes
  int nxn = grd->nxn;
  int nyn = grd->nyn;
  int nzn = grd->nzn;

  // get number of cells
  int nxc = grd->nxc;
  int nyc = grd->nyc;
  int nzc = grd->nzc;

  // allocate temporary arrays
  FPfield ***vectX = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***vectY = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***vectZ = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***imageX = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***imageY = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***imageZ = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***Dx = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***Dy = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***Dz = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***tempX = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***tempY = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***tempZ = newArr3(FPfield, nxn, nyn, nzn);

  FPfield ***divC = newArr3(FPfield, nxc, nyc, nzc);

// set to zero
#pragma omp parallel for
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nzn; k++) {
        vectX[i][j][k] = 0.0;
        vectY[i][j][k] = 0.0;
        vectZ[i][j][k] = 0.0;
        imageX[i][j][k] = 0.0;
        imageY[i][j][k] = 0.0;
        imageZ[i][j][k] = 0.0;
        Dx[i][j][k] = 0.0;
        Dy[i][j][k] = 0.0;
        Dz[i][j][k] = 0.0;
        tempX[i][j][k] = 0.0;
        tempY[i][j][k] = 0.0;
        tempZ[i][j][k] = 0.0;
      }

  // move from krylov space to physical space
  solver2phys(vectX, vectY, vectZ, vector, nxn, nyn, nzn);

  // here we need to impose BC on E: before Laplacian

  ////////
  ////
  //// This part needs to be fixeds. Put PEC
  ///  This puts zero also on the electric field normal
  ///
  //////

  applyBCscalarFieldNzero(vectX, grd, param);
  applyBCscalarFieldNzero(vectY, grd, param);
  applyBCscalarFieldNzero(vectZ, grd, param);

  ////////
  ////
  ////
  ///
  //////

  lapN2N_V(imageX, vectX, grd);
  lapN2N_V(imageY, vectY, grd);
  lapN2N_V(imageZ, vectZ, grd);
  neg(imageX, nxn, nyn, nzn);
  neg(imageY, nxn, nyn, nzn);
  neg(imageZ, nxn, nyn, nzn);

  // grad(div(mu dot E(n + theta)) mu dot E(n + theta) = D
  for (int is = 0; is < param->ns; is++) {
    beta = .5 * param->qom[is] * param->dt / param->c;
#pragma omp parallel for private(omcx, omcy, omcz, edotb, denom)
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
#pragma clang loop vectorize(enable)
        for (int k = 0; k < nzn; k++) {
          omcx = beta * field->Bxn[i][j][k];
          omcy = beta * field->Byn[i][j][k];
          omcz = beta * field->Bzn[i][j][k];
          edotb = vectX[i][j][k] * omcx + vectY[i][j][k] * omcy +
                  vectZ[i][j][k] * omcz;
          denom = param->fourpi / 2 * delt * param->dt / param->c *
                  param->qom[is] * ids[is].rhon[i][j][k] /
                  (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
          Dx[i][j][k] +=
              (vectX[i][j][k] +
               (vectY[i][j][k] * omcz - vectZ[i][j][k] * omcy + edotb * omcx)) *
              denom;
          Dy[i][j][k] +=
              (vectY[i][j][k] +
               (vectZ[i][j][k] * omcx - vectX[i][j][k] * omcz + edotb * omcy)) *
              denom;
          Dz[i][j][k] +=
              (vectZ[i][j][k] +
               (vectX[i][j][k] * omcy - vectY[i][j][k] * omcx + edotb * omcz)) *
              denom;
        }
  }

  // apply boundary condition to Dx, Dy and Dz
  applyBCscalarFieldNzero(Dx, grd, param);
  applyBCscalarFieldNzero(Dy, grd, param);
  applyBCscalarFieldNzero(Dz, grd, param);

  divN2C(divC, Dx, Dy, Dz, grd);

  // communicate you should put BC

  applyBCscalarFieldC(divC, grd, param);

  gradC2N(tempX, tempY, tempZ, divC, grd);

  // -lap(E(n +theta)) - grad(div(mu dot E(n + theta))
  sub(imageX, tempX, nxn, nyn, nzn);
  sub(imageY, tempY, nxn, nyn, nzn);
  sub(imageZ, tempZ, nxn, nyn, nzn);

  // scale delt*delt
  scale(imageX, delt * delt, nxn, nyn, nzn);
  scale(imageY, delt * delt, nxn, nyn, nzn);
  scale(imageZ, delt * delt, nxn, nyn, nzn);

  // -lap(E(n +theta)) - grad(div(mu dot E(n + theta)) + eps dot E(n + theta)
  sum(imageX, Dx, nxn, nyn, nzn);
  sum(imageY, Dy, nxn, nyn, nzn);
  sum(imageZ, Dz, nxn, nyn, nzn);
  sum(imageX, vectX, nxn, nyn, nzn);
  sum(imageY, vectY, nxn, nyn, nzn);
  sum(imageZ, vectZ, nxn, nyn, nzn);

  // move from physical space to krylov space
  phys2solver(im, imageX, imageY, imageZ, nxn, nyn, nzn);

  // deallocate
  delArr3(vectX, nxn, nyn);
  delArr3(vectY, nxn, nyn);
  delArr3(vectZ, nxn, nyn);
  delArr3(imageX, nxn, nyn);
  delArr3(imageY, nxn, nyn);
  delArr3(imageZ, nxn, nyn);
  delArr3(Dx, nxn, nyn);
  delArr3(Dy, nxn, nyn);
  delArr3(Dz, nxn, nyn);
  delArr3(tempX, nxn, nyn);
  delArr3(tempY, nxn, nyn);
  delArr3(tempZ, nxn, nyn);

  delArr3(divC, nxc, nyc);
}

/** calculate the electric field using second order curl-curl formulation of
 * Maxwell equations */
void calculateE(grid *grd, EMfield_aux *field_aux, EMfield *field,
                interpDens_aux *id_aux, interpDensSpecies *ids,
                parameters *param) {
  // get frid points
  int nxn = grd->nxn;
  int nyn = grd->nyn;
  int nzn = grd->nzn;

  int nxc = grd->nxc;
  int nyc = grd->nyc;
  int nzc = grd->nzc;

  // krylov vectors
  FPfield *xkrylov =
      new FPfield[3 * (nxn - 2) * (nyn - 2) * (nzn - 2)];  // 3 E components
  FPfield *bkrylov = new FPfield[3 * (nxn - 2) * (nyn - 2) *
                                 (nzn - 2)];  // 3 components: source

  eqValue(0.0, xkrylov, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2));
  eqValue(0.0, bkrylov, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2));

  std::cout << "*** MAXWELL SOLVER ***" << std::endl;

  // form the source term
  MaxwellSource(bkrylov, grd, field, field_aux, id_aux, param);

  // prepare xkrylov solver
  phys2solver(xkrylov, field->Ex, field->Ey, field->Ez, nxn, nyn, nzn);

  // call the GMRes
  GMRes(&MaxwellImage, xkrylov, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2), bkrylov,
        20, 200, param->GMREStol, field, ids, grd, param);

  // move from krylov space to physical space
  solver2phys(field_aux->Exth, field_aux->Eyth, field_aux->Ezth, xkrylov, nxn,
              nyn, nzn);

  // add E from Eth
  addscale(1 / param->th, -(1.0 - param->th) / param->th, field->Ex,
           field_aux->Exth, nxn, nyn, nzn);
  addscale(1 / param->th, -(1.0 - param->th) / param->th, field->Ey,
           field_aux->Eyth, nxn, nyn, nzn);
  addscale(1 / param->th, -(1.0 - param->th) / param->th, field->Ez,
           field_aux->Ezth, nxn, nyn, nzn);

  // Smooth: You might have special smoothing for E (imposing E on the BC)
  // smoothInterpScalarN(field_aux->Exth, grd, param);
  // smoothInterpScalarN(field_aux->Eyth, grd, param);
  // smoothInterpScalarN(field_aux->Ezth, grd, param);

  // smoothInterpScalarN(field->Ex, grd, param);
  // smoothInterpScalarN(field->Ey, grd, param);
  // smoothInterpScalarN(field->Ez, grd, param);

  // deallocate temporary arrays
  delete[] xkrylov;
  delete[] bkrylov;
}

// calculate the magnetic field from Faraday's law
void calculateB(grid *grd, EMfield_aux *field_aux, EMfield *field,
                parameters *param) {
  int nxc = grd->nxc;
  int nyc = grd->nyc;
  int nzc = grd->nzc;

  FPfield ***tempXC = newArr3(FPfield, nxc, nyc, nzc);
  FPfield ***tempYC = newArr3(FPfield, nxc, nyc, nzc);
  FPfield ***tempZC = newArr3(FPfield, nxc, nyc, nzc);

  std::cout << "*** B CALCULATION ***" << std::endl;

  // calculate the curl of Eth
  curlN2C(tempXC, tempYC, tempZC, field_aux->Exth, field_aux->Eyth,
          field_aux->Ezth, grd);
  // update the magnetic field
  addscale(-param->c * param->dt, 1, field_aux->Bxc, tempXC, nxc, nyc, nzc);
  addscale(-param->c * param->dt, 1, field_aux->Byc, tempYC, nxc, nyc, nzc);
  addscale(-param->c * param->dt, 1, field_aux->Bzc, tempZC, nxc, nyc, nzc);

  applyBCscalarFieldC(field_aux->Bxc, grd, param);
  applyBCscalarFieldC(field_aux->Byc, grd, param);
  applyBCscalarFieldC(field_aux->Bzc, grd, param);

  // interpolate C2N

  interpC2Nfield(field->Bxn, field_aux->Bxc, grd);
  interpC2Nfield(field->Byn, field_aux->Byc, grd);
  interpC2Nfield(field->Bzn, field_aux->Bzc, grd);

  // BC on By ???

  // deallocate
  delArr3(tempXC, nxc, nyc);
  delArr3(tempYC, nxc, nyc);
  delArr3(tempZC, nxc, nyc);
}

/* Poisson Image */
void PoissonImage(FPfield *image, FPfield *vector, grid *grd) {
  // center cells
  int nxc = grd->nxc;
  int nyc = grd->nyc;
  int nzc = grd->nzc;

  // allocate temporary arrays
  FPfield ***temp = newArr3(FPfield, nxc, nyc, nzc);
  FPfield ***im = newArr3(FPfield, nxc, nyc, nzc);

  // set arrays to zero
  for (int i = 0; i < (nxc - 2) * (nyc - 2) * (nzc - 2); i++) image[i] = 0.0;
  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++) {
        temp[i][j][k] = 0.0;
        im[i][j][k] = 0.0;
      }
  solver2phys(temp, vector, nxc, nyc, nzc);
  // the BC for this Laplacian are zero on th eboundary?
  lapC2C(im, temp, grd);
  // move from physical space to krylov space
  phys2solver(image, im, nxc, nyc, nzc);

  // deallocate temporary array and objects
  delArr3(temp, nxc, nyc);
  delArr3(im, nxc, nyc);
}

/** calculate Poisson Correction */
void divergenceCleaning(grid *grd, EMfield_aux *field_aux, EMfield *field,
                        interpDensNet *idn, parameters *param) {
  // get the number of cells
  int nxc = grd->nxc;
  int nyc = grd->nyc;
  int nzc = grd->nzc;
  //  get the number of nodes
  int nxn = grd->nxn;
  int nyn = grd->nyn;
  int nzn = grd->nzn;

  // temporary array for div(E)
  FPfield ***divE = newArr3(FPfield, nxc, nyc, nzc);
  FPfield ***gradPHIX = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***gradPHIY = newArr3(FPfield, nxn, nyn, nzn);
  FPfield ***gradPHIZ = newArr3(FPfield, nxn, nyn, nzn);

  // 1D vectors for solver
  FPfield *xkrylovPoisson = new FPfield[(nxc - 2) * (nyc - 2) * (nzc - 2)];
  FPinterp *bkrylovPoisson = new FPinterp[(nxc - 2) * (nyc - 2) * (nzc - 2)];
  // set to zero xkrylov and bkrylov
  for (int i = 0; i < ((nxc - 2) * (nyc - 2) * (nzc - 2)); i++) {
    xkrylovPoisson[i] = 0.0;
    bkrylovPoisson[i] = 0.0;
  }
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) {
        gradPHIX[i][j][k] = 0.0;
        gradPHIY[i][j][k] = 0.0;
        gradPHIZ[i][j][k] = 0.0;
      }

  std::cout << "*** DIVERGENCE CLEANING ***" << std::endl;
  divN2C(divE, field->Ex, field->Ey, field->Ez, grd);

  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++) {
        divE[i][j][k] = divE[i][j][k] - param->fourpi * idn->rhoc[i][j][k];
        field_aux->Phi[i][j][k] = 0.0;
      }

  // phys to solver
  phys2solver(bkrylovPoisson, divE, nxc, nyc, nzc);
  // call CG for solving Poisson
  if (!CG(xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2), bkrylovPoisson,
          3000, param->CGtol, &PoissonImage, grd))
    std::cout << "*ERROR - CG not Converged" << std::endl;
  solver2phys(field_aux->Phi, xkrylovPoisson, nxc, nyc, nzc);

  // This has Newmann. If commented has zero in ghost cells
  applyBCscalarFieldC(field_aux->Phi, grd, param);
  gradC2N(gradPHIX, gradPHIY, gradPHIZ, field_aux->Phi, grd);

  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) {
        field->Ex[i][j][k] -= gradPHIX[i][j][k];
        field->Ey[i][j][k] -= gradPHIY[i][j][k];
        field->Ez[i][j][k] -= gradPHIZ[i][j][k];
      }

  // deallocate vectors
  delete[] xkrylovPoisson;
  delete[] bkrylovPoisson;

  // deallocate temporary array
  delArr3(divE, nxc, nyc);
  delArr3(gradPHIX, nxn, nyn);
  delArr3(gradPHIY, nxn, nyn);
  delArr3(gradPHIZ, nxn, nyn);
}

/** Calculate hat densities: Jh and rhoh*/
void calculateHatDensities(struct interpDens_aux *id_aux,
                           struct interpDensNet *idn,
                           struct interpDensSpecies *ids, struct EMfield *field,
                           struct grid *grd, struct parameters *param) {
  // parameters
  FPfield beta, edotb, omcx, omcy, omcz, denom;

  // nodes
  int nxn = grd->nxn;
  int nyn = grd->nyn;
  int nzn = grd->nzn;

  // centers
  int nxc = grd->nxc;
  int nyc = grd->nyc;
  int nzc = grd->nzc;

  // allocate temporary ararys
  // center
  FPinterp ***tempXC = newArr3(FPinterp, nxc, nyc, nzc);
  FPinterp ***tempYC = newArr3(FPinterp, nxc, nyc, nzc);
  FPinterp ***tempZC = newArr3(FPinterp, nxc, nyc, nzc);
  // nodes
  FPinterp ***tempXN = newArr3(FPinterp, nxn, nyn, nzn);
  FPinterp ***tempYN = newArr3(FPinterp, nxn, nyn, nzn);
  FPinterp ***tempZN = newArr3(FPinterp, nxn, nyn, nzn);

  // Set J hat to zero: Important to set to zero because then it accumulate in
  // PIdot
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) {
        id_aux->Jxh[i][j][k] = 0.0;
        id_aux->Jyh[i][j][k] = 0.0;
        id_aux->Jzh[i][j][k] = 0.0;
      }
  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++) id_aux->rhoh[i][j][k] = 0.0;

  // smoothing of rhoc: BC have applied before returning
  // smoothInterpScalarC(idn->rhoc,grd,param);
  // apply boundary to rhoh
  applyBCscalarDensC(id_aux->rhoh, grd,
                     param);  // set BC on ghost cells before interp

  for (int is = 0; is < param->ns; is++) {
    divSymmTensorN2C(tempXC, tempYC, tempZC, ids[is].pxx, ids[is].pxy,
                     ids[is].pxz, ids[is].pyy, ids[is].pyz, ids[is].pzz, grd);

    // scale the pressure tensor
    scale(tempXC, -param->dt / 2.0, nxc, nyc, nzc);
    scale(tempYC, -param->dt / 2.0, nxc, nyc, nzc);
    scale(tempZC, -param->dt / 2.0, nxc, nyc, nzc);

    // apply BC to centers before interpolation: this is not needed
    applyBCscalarDensC(tempXC, grd,
                       param);  // set BC on ghost cells before interp
    applyBCscalarDensC(tempYC, grd,
                       param);  // set BC on ghost cells before interp
    applyBCscalarDensC(tempZC, grd,
                       param);  // set BC on ghost cells before interp

    // interpolation
    interpC2Ninterp(tempXN, tempXC, grd);
    interpC2Ninterp(tempYN, tempYC, grd);
    interpC2Ninterp(tempZN, tempZC, grd);

    // sum(tempXN, Jxs, nxn, nyn, nzn, is); sum(tempYN, Jys, nxn, nyn, nzn, is);
    // sum(tempZN, Jzs, nxn, nyn, nzn, is);
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          tempXN[i][j][k] += ids[is].Jx[i][j][k];
          tempYN[i][j][k] += ids[is].Jy[i][j][k];
          tempZN[i][j][k] += ids[is].Jz[i][j][k];
        }

    // PIDOT //PIdot(Jxh, Jyh, Jzh, tempXN, tempYN, tempZN, is, grid);
    beta = .5 * param->qom[is] * param->dt / param->c;
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          omcx = beta * field->Bxn[i][j][k];
          omcy = beta * field->Byn[i][j][k];
          omcz = beta * field->Bzn[i][j][k];
          edotb = tempXN[i][j][k] * omcx + tempYN[i][j][k] * omcy +
                  tempZN[i][j][k] * omcz;
          denom = 1 / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
          id_aux->Jxh[i][j][k] +=
              (tempXN[i][j][k] + (tempYN[i][j][k] * omcz -
                                  tempZN[i][j][k] * omcy + edotb * omcx)) *
              denom;
          id_aux->Jyh[i][j][k] +=
              (tempYN[i][j][k] + (tempZN[i][j][k] * omcx -
                                  tempXN[i][j][k] * omcz + edotb * omcy)) *
              denom;
          id_aux->Jzh[i][j][k] +=
              (tempZN[i][j][k] + (tempXN[i][j][k] * omcy -
                                  tempYN[i][j][k] * omcx + edotb * omcz)) *
              denom;
        }
  }

  // smooth J hat
  // smoothInterpScalarN(id_aux->Jxh,grd,param);
  // smoothInterpScalarN(id_aux->Jyh,grd,param);
  // smoothInterpScalarN(id_aux->Jzh,grd,param);

  // calculate rho hat = rho - (dt*theta)div(jhat)
  // set tempXC to zero
  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++) tempXC[i][j][k] = 0.0;

  // in principle dont need this
  applyBCscalarDensN(id_aux->Jxh, grd, param);
  applyBCscalarDensN(id_aux->Jyh, grd, param);
  applyBCscalarDensN(id_aux->Jzh, grd, param);

  divN2C(tempXC, id_aux->Jxh, id_aux->Jyh, id_aux->Jzh, grd);

  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
        id_aux->rhoh[i][j][k] =
            idn->rhoc[i][j][k] - param->dt * param->th * tempXC[i][j][k];

  // apply boundary to rhoh
  applyBCscalarDensC(id_aux->rhoh, grd,
                     param);  // set BC on ghost cells before interp

  // deallocate
  delArr3(tempXC, nxc, nyc);
  delArr3(tempYC, nxc, nyc);
  delArr3(tempZC, nxc, nyc);
  delArr3(tempXN, nxn, nyn);
  delArr3(tempYN, nxn, nyn);
  delArr3(tempZN, nxn, nyn);
}
