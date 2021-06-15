#ifndef TRACKING_H
#define TRACKING_H

#include "Parameters.h"
#include "Particles.h"
#include "PrecisionTypes.h"

#include <vector>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <iostream>

void find_and_toggle_track_particles(struct parameters* param,
                                        struct particles* part);

#endif