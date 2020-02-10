#include "Particles.h"
#include <math.h>

/** allocate particle arrays */
void particle_allocate(struct parameters *param, struct particles *part,
                       int is) {
  // set species ID
  part->species_ID = is;
  // number of particles
  part->nop = param->np[is];
  // maximum number of particles
  part->npmax = param->npMax[is];

  // choose a different number of mover iterations for ions and electrons
  part->NiterMover = param->NiterMover;
  part->n_sub_cycles = param->n_sub_cycles;

  // particles per cell
  part->npcelx = param->npcelx[is];
  part->npcely = param->npcely[is];
  part->npcelz = param->npcelz[is];
  part->npcel = part->npcelx * part->npcely * part->npcelz;

  // cast it to required precision
  part->qom = (FPpart)param->qom[is];

  long npmax = part->npmax;

  // initialize drift and thermal velocities
  // drift
  part->u0 = (FPpart)param->u0[is];
  part->v0 = (FPpart)param->v0[is];
  part->w0 = (FPpart)param->w0[is];
  // thermal
  part->uth = (FPpart)param->uth[is];
  part->vth = (FPpart)param->vth[is];
  part->wth = (FPpart)param->wth[is];

  //////////////////////////////
  /// ALLOCATION PARTICLE ARRAYS
  //////////////////////////////
  part->x = new FPpart[npmax];
  part->y = new FPpart[npmax];
  part->z = new FPpart[npmax];
  // allocate velocity
  part->u = new FPpart[npmax];
  part->v = new FPpart[npmax];
  part->w = new FPpart[npmax];
  // allocate charge = q * statistical weight
  part->q = new FPinterp[npmax];

  part->track_particle = new bool[npmax];
}
/** deallocate */
void particle_deallocate(struct particles *part) {
  // deallocate particle variables
  delete[] part->x;
  delete[] part->y;
  delete[] part->z;
  delete[] part->u;
  delete[] part->v;
  delete[] part->w;
  delete[] part->q;
}

/** allocate particle arrays */
void particle_aux_allocate(struct particles *part,
                           struct particles_aux *part_aux, int is) {
  // set species ID
  part_aux->species_ID = is;
  // number of particles
  part_aux->nop = part->nop;
  // maximum number of particles
  part_aux->npmax = part->npmax;

  long npmax = part->npmax;

  // allocate densities brought by each particle
  part_aux->rho_p = new FPpart[part->npmax][2][2][2];
  part_aux->Jx = new FPpart[part->npmax][2][2][2];
  part_aux->Jy = new FPpart[part->npmax][2][2][2];
  part_aux->Jz = new FPpart[part->npmax][2][2][2];
  part_aux->pxx = new FPpart[part->npmax][2][2][2];
  part_aux->pxy = new FPpart[part->npmax][2][2][2];
  part_aux->pxz = new FPpart[part->npmax][2][2][2];
  part_aux->pyy = new FPpart[part->npmax][2][2][2];
  part_aux->pyz = new FPpart[part->npmax][2][2][2];
  part_aux->pzz = new FPpart[part->npmax][2][2][2];

  // cell index
  part_aux->ix_p = new int[part->npmax];
  part_aux->iy_p = new int[part->npmax];
  part_aux->iz_p = new int[part->npmax];
}

void particle_aux_deallocate(struct particles_aux *part_aux) {
  // deallocate auxiliary particle variables needed for particle interpolation
  delete[] part_aux->rho_p;
  delete[] part_aux->Jx;
  delete[] part_aux->Jy;
  delete[] part_aux->Jz;
  delete[] part_aux->pxx;
  delete[] part_aux->pxy;
  delete[] part_aux->pxz;
  delete[] part_aux->pyy;
  delete[] part_aux->pyz;
  delete[] part_aux->pzz;
}
