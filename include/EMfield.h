#ifndef EMFIELD_H
#define EMFIELD_H

#include "Grid.h"
#include "PrecisionTypes.h"

/** structure with field information */
struct EMfield {
  // field arrays: 4D arrays

  /* Electric field defined on nodes: last index is component */
  FPfield ***Ex;
  FPfield *Ex_flat;
  FPfield ***Ey;
  FPfield *Ey_flat;
  FPfield ***Ez;
  FPfield *Ez_flat;
  /* Magnetic field defined on nodes: last index is component */
  FPfield ***Bxn;
  FPfield *Bxn_flat;
  FPfield ***Byn;
  FPfield *Byn_flat;
  FPfield ***Bzn;
  FPfield *Bzn_flat;
};

/** structure with auxiliary field quantities like potentials or quantities
 * defined at centers  */
struct EMfield_aux {
  /* Electrostatic potential defined on central points*/
  FPfield ***Phi;
  FPfield *Phi_flat;

  /* Electric field at time theta */
  FPfield ***Exth;
  FPfield *Exth_flat;

  FPfield ***Eyth;
  FPfield *Eyth_flat;

  FPfield ***Ezth;
  FPfield *Ezth_flat;

  /* Magnetic field defined on nodes: last index is component - Centers */
  FPfield ***Bxc;
  FPfield *Bxc_flat;
  FPfield ***Byc;
  FPfield *Byc_flat;
  FPfield ***Bzc;
  FPfield *Bzc_flat;
};

/** allocate electric and magnetic field */
void field_allocate(struct grid *grd, struct EMfield *field);

/** deallocate electric and magnetic field */
void field_deallocate(struct grid *grd, struct EMfield *field);

/** allocate electric and magnetic field */
void field_aux_allocate(struct grid *grd, struct EMfield_aux *field_aux);

/** deallocate */
void field_aux_deallocate(struct grid *grd, struct EMfield_aux *field_aux);

#endif
