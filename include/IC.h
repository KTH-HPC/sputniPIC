#ifndef IC_H
#define IC_H

#include <math.h>

#include "EMfield.h"
#include "Grid.h"
#include "InterpDensSpecies.h"
#include "Parameters.h"
#include "Particles.h"

/** initialize for magnetic reconnection probelm with Harris current sheet */
void initGEM(struct parameters *param, struct grid *grd, struct EMfield *field,
             struct EMfield_aux *field_aux, struct particles *part,
             struct interpDensSpecies *ids);

/** initialize uniform electric and magnetic field */
void initUniform(struct parameters *param, struct grid *grd,
                 struct EMfield *field, struct EMfield_aux *field_aux,
                 struct particles *part, struct interpDensSpecies *ids);

void read_ic_data(struct EMfield *field, struct EMfield_aux *field_aux, struct grid *grd, struct interpDensSpecies *ids, struct particles *part, struct parameters *param);
void save_ic_data(struct EMfield *field, struct EMfield_aux *field_aux, struct grid *grd, struct interpDensSpecies *ids, struct particles *part, struct parameters *param);

#endif
