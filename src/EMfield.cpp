#include "EMfield.h"

/** allocate electric and magnetic field */
void field_allocate(struct grid *grd, struct EMfield *field) {
  // E on nodes
  field->Ex = newArr3<FPfield>(&field->Ex_flat, grd->nxn, grd->nyn, grd->nzn);
  field->Ey = newArr3<FPfield>(&field->Ey_flat, grd->nxn, grd->nyn, grd->nzn);
  field->Ez = newArr3<FPfield>(&field->Ez_flat, grd->nxn, grd->nyn, grd->nzn);
  // B on nodes
  field->Bxn = newArr3<FPfield>(&field->Bxn_flat, grd->nxn, grd->nyn, grd->nzn);
  field->Byn = newArr3<FPfield>(&field->Byn_flat, grd->nxn, grd->nyn, grd->nzn);
  field->Bzn = newArr3<FPfield>(&field->Bzn_flat, grd->nxn, grd->nyn, grd->nzn);
}

/** deallocate electric and magnetic field */
void field_deallocate(struct grid *grd, struct EMfield *field) {
  // E deallocate 3D arrays
  delArr3(field->Ex, grd->nxn, grd->nyn);
  delArr3(field->Ey, grd->nxn, grd->nyn);
  delArr3(field->Ez, grd->nxn, grd->nyn);
  // B deallocate 3D arrays
  delArr3(field->Bxn, grd->nxn, grd->nyn);
  delArr3(field->Byn, grd->nxn, grd->nyn);
  delArr3(field->Bzn, grd->nxn, grd->nyn);
}

/** allocate electric and magnetic field */
void field_aux_allocate(struct grid *grd, struct EMfield_aux *field_aux) {
  // Electrostatic potential
  field_aux->Phi =
      newArr3<FPfield>(&field_aux->Phi_flat, grd->nxc, grd->nyc, grd->nzc);

  // allocate 3D arrays
  field_aux->Exth =
      newArr3<FPfield>(&field_aux->Exth_flat, grd->nxn, grd->nyn, grd->nzn);
  field_aux->Eyth =
      newArr3<FPfield>(&field_aux->Eyth_flat, grd->nxn, grd->nyn, grd->nzn);
  field_aux->Ezth =
      newArr3<FPfield>(&field_aux->Ezth_flat, grd->nxn, grd->nyn, grd->nzn);
  // B on centers
  field_aux->Bxc =
      newArr3<FPfield>(&field_aux->Bxc_flat, grd->nxc, grd->nyc, grd->nzc);
  field_aux->Byc =
      newArr3<FPfield>(&field_aux->Byc_flat, grd->nxc, grd->nyc, grd->nzc);
  field_aux->Bzc =
      newArr3<FPfield>(&field_aux->Bzc_flat, grd->nxc, grd->nyc, grd->nzc);
}

/** deallocate */
void field_aux_deallocate(struct grid *grd, struct EMfield_aux *field_aux) {
  // Eth
  delArr3(field_aux->Exth, grd->nxn, grd->nyn);
  delArr3(field_aux->Eyth, grd->nxn, grd->nyn);
  delArr3(field_aux->Ezth, grd->nxn, grd->nyn);
  // Bc
  delArr3(field_aux->Bxc, grd->nxc, grd->nyc);
  delArr3(field_aux->Byc, grd->nxc, grd->nyc);
  delArr3(field_aux->Bzc, grd->nxc, grd->nyc);
  // Phi
  delArr3(field_aux->Phi, grd->nxc, grd->nyc);
}
