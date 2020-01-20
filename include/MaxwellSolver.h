#ifndef MAXWELLSOLVER
#define MAXWELLSOLVER

#include "Basic.h"
#include "EMfield.h"
#include "InterpDensNet.h"
#include "InterpDensSpecies.h"
#include "InterpDens_aux.h"
#include "PrecisionTypes.h"
#include "Solvers.h"

void MaxwellImage(FPfield *im, FPfield *vector, EMfield *field,
                  interpDensSpecies *ids, grid *grd, parameters *param);

void MaxwellSource(FPfield *bkrylov, grid *grd, EMfield *field,
                   EMfield_aux *field_aux, interpDens_aux *id_aux,
                   parameters *param);

/** calculate the electric field using second order curl-curl formulation of
 * Maxwell equations */
void calculateE(grid *grd, EMfield_aux *field_aux, EMfield *field,
                interpDens_aux *id_aux, interpDensSpecies *ids,
                parameters *param);

// calculate the magnetic field from Faraday's law
void calculateB(grid *grd, EMfield_aux *field_aux, EMfield *field,
                parameters *param);

/* Poisson Image */
void PoissonImage(FPfield *image, FPfield *vector, grid *grd);

/** calculate Poisson Correction */
void divergenceCleaning(grid *grd, EMfield_aux *field_aux, EMfield *field,
                        interpDensNet *idn, parameters *param);

/** Calculate hat densities: Jh and rhoh*/
void calculateHatDensities(struct interpDens_aux *id_aux,
                           struct interpDensNet *idn,
                           struct interpDensSpecies *ids, struct EMfield *field,
                           struct grid *grd, struct parameters *param);

#endif
