#include "Tracking.h"

void find_and_toggle_track_particles(struct parameters* param,
                                     struct particles* part) {
  FPpart center_x = param->Lx / 2.0;
  FPpart center_y = param->Ly / 2.0;
  FPpart center_z = param->Lz / 2.0;
  FPpart min_x = center_x - (param->tracking_Lx / 2.0);
  FPpart max_x = center_x + (param->tracking_Lx / 2.0);
  FPpart min_y = center_y - (param->tracking_Ly / 2.0);
  FPpart max_y = center_y + (param->tracking_Ly / 2.0);
  FPpart min_z = center_z - (param->tracking_Lz / 2.0);
  FPpart max_z = center_z + (param->tracking_Lz / 2.0);

  std::vector<size_t> particles_in_center;
  // find all particles that are in the center area.
  for (size_t i = 0; i < part->npmax; i++) {
    if (part->x[i] > min_x && part->x[i] < max_x && part->y[i] > min_y &&
        part->y[i] < max_y && part->z[i] > min_z && part->z[i] < max_z) {
      particles_in_center.push_back(i);
    }
  }
  particles_in_center.shrink_to_fit();

  if (particles_in_center.size() < param->n_tracked_particles) {
    std::cerr << "Fewer particles than requested were sampled for tracking!"
              << std::endl;
  }

  // shuffle list of tracked particles.
  auto rng = std::default_random_engine{};
  std::shuffle(particles_in_center.begin(), particles_in_center.end(), rng);

  // select n first elements.
  particles_in_center.resize(param->n_tracked_particles);
  for (size_t& i : particles_in_center) {
    part->track_particle[i] = true;
  }
}
