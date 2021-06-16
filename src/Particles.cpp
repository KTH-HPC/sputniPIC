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
  std::fill_n(part->track_particle, npmax, 0);
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

inline long bin_index(long x, long y, long z, long ny, long nz) {
  return z + y * nz + x * nz * ny;
}

/** sort particles with regards to grid layout **/
void particle_sort(struct parameters *params, struct particles *particles,
                   struct grid *grd) {
  /**
   * This function can be done in two ways.
   * The first is to simply loop through the sorting zones one by one
   * and go through all particles for each zone and adding them to the
   * results as they are checked. This approach is simple/low memory but slow.
   * The other approach is to setup some sort of bins (map of std::vector?)
   * where sorted particle indexes are stored which means that the sorting
   * only needs one pass over the particles. The downside is potentially worse
   * worst case performance plus more overhead.
   */

#ifdef NAIVE_SORT
  // Loop over every bin and add all particles matching to the new arrays.
  // Add the particles not found when looping over the bins (edge cases).

  std::vector<bool> is_particle_sorted(particles->nop, false);

  FPpart *x = new FPpart[particles->npmax];
  FPpart *y = new FPpart[particles->npmax];
  FPpart *z = new FPpart[particles->npmax];
  FPpart *u = new FPpart[particles->npmax];
  FPpart *v = new FPpart[particles->npmax];
  FPpart *w = new FPpart[particles->npmax];
  FPinterp *q = new FPpart[particles->npmax];
  long p_counter = 0;

  // Calculate and setup number of sorting zones.
  for (long x_zone = 0; x_zone < grd->nxc; x_zone += params->sort_cps) {
    for (long y_zone = 0; y_zone < grd->nyc; y_zone += params->sort_cps) {
      for (long z_zone = 0; z_zone < grd->nzc; z_zone += params->sort_cps) {
        for (long i = 0; i < particles->nop; i++) {
          if (!is_particle_sorted[i]) {
            if (particles->x[i] >= x_zone * grd->dx &&
                particles->x[i] < (x_zone + params->sort_cps + 1) * grd->dx &&
                particles->y[i] >= y_zone * grd->dy &&
                particles->y[i] < (y_zone + params->sort_cps + 1) * grd->dy &&
                particles->z[i] >= z_zone * grd->dz &&
                particles->z[i] < (z_zone + params->sort_cps + 1) * grd->dz) {
              is_particle_sorted[i] = true;
              x[p_counter] = particles->x[p_counter];
              y[p_counter] = particles->y[p_counter];
              z[p_counter] = particles->z[p_counter];
              u[p_counter] = particles->u[p_counter];
              v[p_counter] = particles->v[p_counter];
              w[p_counter] = particles->w[p_counter];
              q[p_counter] = particles->q[p_counter];
              p_counter++;
            }
          }
        }
      }
    }
  }

  // One last loop to catch the particles missed earlier.
  for (long i = 0; i < particles->nop; i++) {
    if (!is_particle_sorted[i]) {
      x[p_counter] = particles->x[p_counter];
      y[p_counter] = particles->y[p_counter];
      z[p_counter] = particles->z[p_counter];
      u[p_counter] = particles->u[p_counter];
      v[p_counter] = particles->v[p_counter];
      w[p_counter] = particles->w[p_counter];
      q[p_counter] = particles->q[p_counter];
      p_counter++;
    }
  }

  std::memcpy(particles->x, x, particles->npmax * sizeof(x[0]));
  std::memcpy(particles->y, y, particles->npmax * sizeof(y[0]));
  std::memcpy(particles->z, z, particles->npmax * sizeof(z[0]));
  std::memcpy(particles->u, u, particles->npmax * sizeof(u[0]));
  std::memcpy(particles->v, v, particles->npmax * sizeof(v[0]));
  std::memcpy(particles->w, w, particles->npmax * sizeof(w[0]));
  std::memcpy(particles->q, q, particles->npmax * sizeof(q[0]));

  delete[] x;
  delete[] y;
  delete[] z;
  delete[] u;
  delete[] v;
  delete[] w;
  delete[] q;
#else
  // Loop over every particle and add them to the corresponding bin.
  // For each bin, add to the new particle arrays in order.
  long n_binx = 1 + ((grd->nxc - 1) / params->sort_cps);
  long n_biny = 1 + ((grd->nyc - 1) / params->sort_cps);
  long n_binz = 1 + ((grd->nzc - 1) / params->sort_cps);
  FPpart dist_binx = grd->Lx / (FPpart)n_binx;
  FPpart dist_biny = grd->Ly / (FPpart)n_biny;
  FPpart dist_binz = grd->Lz / (FPpart)n_binz;

  long bins_tot = n_binx * n_biny * n_binz;

  // Create a vector containing vectors that represent each bin where
  // particle indexes can be stored. Indexing is done z -> y -> x.
  // Add extra bin for "lost" particles.
  std::vector<std::vector<long> > bins(bins_tot + 1, std::vector<long>());

  // Reserve space ahead of usage.
  for (auto bin : bins) {
    bin.reserve(particles->npmax / bins_tot);
  }

  for (long i = 0; i < particles->npmax; i++) {
    // Calculate which bin the particle should be in.
    long pos_x = floor(particles->x[i] / dist_binx);
    long pos_y = floor(particles->y[i] / dist_biny);
    long pos_z = floor(particles->z[i] / dist_binz);

    if (pos_x < 0 || pos_x >= n_binx || pos_y < 0 || pos_y >= n_biny ||
        pos_z < 0 || pos_z >= n_binz) {
      // Particle is not inside the bin zone.
      bins.back().push_back(i);
    } else {
      bins[bin_index(pos_x, pos_y, pos_z, n_biny, n_binz)].push_back(i);
    }
  }

  // Transfer the particle information to the new arrays.
  FPpart *x = new FPpart[particles->npmax];
  FPpart *y = new FPpart[particles->npmax];
  FPpart *z = new FPpart[particles->npmax];
  FPpart *u = new FPpart[particles->npmax];
  FPpart *v = new FPpart[particles->npmax];
  FPpart *w = new FPpart[particles->npmax];
  FPinterp *q = new FPpart[particles->npmax];

  long p_counter = 0;
  for (auto bin : bins) {
    for (auto part_pos : bin) {
      x[p_counter] = particles->x[part_pos];
      y[p_counter] = particles->y[part_pos];
      z[p_counter] = particles->z[part_pos];
      u[p_counter] = particles->u[part_pos];
      v[p_counter] = particles->v[part_pos];
      w[p_counter] = particles->w[part_pos];
      q[p_counter] = particles->q[part_pos];
      p_counter++;
    }
  }

  std::memcpy(particles->x, x, particles->npmax * sizeof(x[0]));
  std::memcpy(particles->y, y, particles->npmax * sizeof(y[0]));
  std::memcpy(particles->z, z, particles->npmax * sizeof(z[0]));
  std::memcpy(particles->u, u, particles->npmax * sizeof(u[0]));
  std::memcpy(particles->v, v, particles->npmax * sizeof(v[0]));
  std::memcpy(particles->w, w, particles->npmax * sizeof(w[0]));
  std::memcpy(particles->q, q, particles->npmax * sizeof(q[0]));

  delete[] x;
  delete[] y;
  delete[] z;
  delete[] u;
  delete[] v;
  delete[] w;
  delete[] q;
#endif
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
