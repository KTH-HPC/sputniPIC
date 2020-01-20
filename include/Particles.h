#ifndef PARTICLES_H
#define PARTICLES_H

#include <math.h>
#include "Parameters.h"
#include "PrecisionTypes.h"

struct particles {
  /** species ID: 0, 1, 2 , ... */
  int species_ID;

  /** maximum number of particles of this species on this domain. used for
   * memory allocation */
  long npmax;
  /** number of particles of this species on this domain */
  long nop;

  /** Electron and ions have different number of iterations: ions moves slower
   * than ions */
  int NiterMover;
  /** number of particle of subcycles in the mover */
  int n_sub_cycles;

  /** number of particles per cell */
  int npcel;
  /** number of particles per cell - X direction */
  int npcelx;
  /** number of particles per cell - Y direction */
  int npcely;
  /** number of particles per cell - Z direction */
  int npcelz;

  /** charge over mass ratio */
  FPpart qom;

  /* drift and thermal velocities for this species */
  FPpart u0, v0, w0;
  FPpart uth, vth, wth;

  /** particle arrays: 1D arrays[npmax] */
  FPpart *x;
  FPpart *y;
  FPpart *z;
  FPpart *u;
  FPpart *v;
  FPpart *w;
  /** q must have precision of interpolated quantities: typically double. Not
   * used in mover */
  FPinterp *q;

  /** tracked particle*/
  bool *track_particle;
};

struct particles_aux {
  /** species ID: 0, 1, 2 , ... */
  int species_ID;

  /** maximum number of particles of this species on this domain. used for
   * memory allocation */
  long npmax;
  /** number of particles of this species on this domain */
  long nop;

  /** densities carried nop,2,2,2*/
  FPpart (*rho_p)[2][2][2];
  FPpart (*Jx)[2][2][2];
  FPpart (*Jy)[2][2][2];
  FPpart (*Jz)[2][2][2];
  FPpart (*pxx)[2][2][2];
  FPpart (*pxy)[2][2][2];
  FPpart (*pxz)[2][2][2];
  FPpart (*pyy)[2][2][2];
  FPpart (*pyz)[2][2][2];
  FPpart (*pzz)[2][2][2];

  /** cell index: ix, iy, iz */
  int *ix_p;
  int *iy_p;
  int *iz_p;
};

/** allocate particle arrays */
void particle_allocate(struct parameters *param, struct particles *part,
                       int is);

/** deallocate */
void particle_deallocate(struct particles *part);

void particle_aux_allocate(struct particles *part,
                           struct particles_aux *part_aux, int is);
void particle_aux_deallocate(struct particles_aux *part_aux);

#endif
