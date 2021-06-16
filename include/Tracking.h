#ifndef TRACKING_H
#define TRACKING_H

#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include "Parameters.h"
#include "Particles.h"
#include "PrecisionTypes.h"

void find_and_toggle_track_particles(struct parameters* param,
                                     struct particles* part);

#endif