#ifndef INTERPDENS_AUX_H
#define INTERPDENS_AUX_H

#include "Smoothing.h"

struct interpDens_aux {
  /** charged densities */
  FPinterp ***rhoh;  // rho hat defined at center cell
  /** J current densities */
  FPinterp ***Jxh;
  FPinterp ***Jyh;
  FPinterp ***Jzh;  // on nodes
};

/** allocated interpolated densities per species */
inline void interp_dens_aux_allocate(struct grid *grd,
                                     struct interpDens_aux *id_aux) {
  // hat - charge density defined on nodes and center cell
  id_aux->rhoh = newArr3(FPinterp, grd->nxc, grd->nyc, grd->nzc);
  // hat - current
  id_aux->Jxh = newArr3(FPinterp, grd->nxn, grd->nyn, grd->nzn);
  id_aux->Jyh = newArr3(FPinterp, grd->nxn, grd->nyn, grd->nzn);
  id_aux->Jzh = newArr3(FPinterp, grd->nxn, grd->nyn, grd->nzn);
}

/** deallocate interpolated densities per species */
inline void interp_dens_aux_deallocate(struct grid *grd,
                                       struct interpDens_aux *id_aux) {
  // hat - charge density
  delArr3(id_aux->rhoh, grd->nxc, grd->nyc);
  // hat - current
  delArr3(id_aux->Jxh, grd->nxn, grd->nyn);
  delArr3(id_aux->Jyh, grd->nxn, grd->nyn);
  delArr3(id_aux->Jzh, grd->nxn, grd->nyn);
}

#endif
