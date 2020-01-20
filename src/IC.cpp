#include "IC.h"
#include <math.h>
#include <iostream>

/** initialize for magnetic reconnection probelm with Harris current sheet */
void initGEM(struct parameters* param, struct grid* grd, struct EMfield* field, struct EMfield_aux* field_aux, struct particles* part, struct interpDensSpecies* ids){
    
    // perturbation localized in X
    double pertX = 0.4;
    double xpert, ypert, exp_pert;
    
    // print settings
    std::cout << "*************************************************" << std::endl;
    std::cout << "**  Initialize GEM Challenge with Pertubation  **" << std::endl;
    std::cout << "*************************************************" << std::endl;
    std::cout << "** B0x = " << param->B0x << std::endl;
    std::cout << "** B0y = " << param->B0y << std::endl;
    std::cout << "** B0z = " << param->B0z << std::endl;
    std::cout << "** Delta (current sheet thickness) = " << param->delta << std::endl;
    for (int is=0; is < param->ns; is++){
        std::cout << "** rho species " << is << " = " << param->rhoINIT[is];
        if (is < 2)
            std::cout << " CURRENT SHEET " << std::endl;
        else
            std::cout << " BACKGROUND " << std::endl;
    }
    std::cout << "*************************************************" << std::endl;
    
    /////////////////////////////////////////////////
    //////   FIELD AND DENSITY
    /////////////////////////////////////////////////
    // Set the electric field, magnetic field + rhos
    for (int i=0; i < grd->nxn; i++)
        for (int j=0; j < grd->nyn; j++)
            for (int k=0; k < grd->nzn; k++){
                // initialize the density for species
                for (int is=0; is < param->ns; is++){
                    if (is < 2) // current sheet
                        ids[is].rhon[i][j][k] = (FPinterp) ((param->rhoINIT[is] / (cosh((grd->YN[i][j][k] - grd->Ly/2)/ param->delta) *cosh((grd->YN[i][j][k] - grd->Ly/2) / param->delta)))) / param->fourpi;
                    else // background
                        ids[is].rhon[i][j][k] = (FPinterp) param->rhoINIT[is]/param->fourpi;
                }
                //std::cout << "OK" << std::endl;
                // electric field
                // electric field
                field->Ex[i][j][k] =  0.0;
                field->Ey[i][j][k] =  0.0;
                field->Ez[i][j][k] =  0.0;
                field_aux->Exth[i][j][k] =  0.0;
                field_aux->Eyth[i][j][k] =  0.0;
                field_aux->Ezth[i][j][k] =  0.0;
                // Magnetic field
                field->Bxn[i][j][k] = param->B0x*tanh((grd->YN[i][j][k] - grd->Ly/2)/param->delta);
                // add the initial GEM perturbation
                // Bxn[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXN(i,j,k)/Lx)*sin(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly  );
                field->Byn[i][j][k] = param->B0y; // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXN(i,j,k)/Lx)*cos(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly);
                // add the initial X perturbation
                xpert = grd->XN[i][j][k]- grd->Lx/2;
                ypert = grd->YN[i][j][k]- grd->Ly/2;
                exp_pert = exp(-(xpert/param->delta)*(xpert/param->delta)-(ypert/param->delta)*(ypert/param->delta));
                field->Bxn[i][j][k] +=(param->B0x*pertX)*exp_pert*(-cos(M_PI*xpert/10.0/param->delta)*cos(M_PI*ypert/10.0/param->delta)*2.0*ypert/param->delta
                                                     -cos(M_PI*xpert/10.0/param->delta)*sin(M_PI*ypert/10.0/param->delta)*M_PI/10.0);
                field->Byn[i][j][k] +=(param->B0x*pertX)*exp_pert*(cos(M_PI*xpert/10.0/param->delta)*cos(M_PI*ypert/10.0/param->delta)*2.0*xpert/param->delta
                                                     +sin(M_PI*xpert/10.0/param->delta)*cos(M_PI*ypert/10.0/param->delta)*M_PI/10.0);
                // guide field
                field->Bzn[i][j][k] = param->B0z;
                
            }
    // calculate B and rho at centers cells: first argument is on center cell
    interpN2Cfield(field_aux->Bxc, field_aux->Byc, field_aux->Bzc, field->Bxn, field->Byn, field->Bzn, grd);
    

    // interpolate densities species from node
    for (int is=0; is < param->ns; is++){
        interpN2Crho(&ids[is], grd);
    }
    
    /////////////////////////////////////////////////
    //////   PARTICLE
    /////////////////////////////////////////////////
    
    double harvest;
    double prob, theta, sign;
    long long counter;
    
    // loop over the species
    for (int is=0; is < param->ns; is++){
        // set particle counter to zero
        counter=0;
        // set the seed for random number generator equal to species id
        srand(is);
        for (int i=1; i< grd->nxc-1;i++)
            for (int j=1; j< grd->nyc-1;j++)
                for (int k=1; k< grd->nzc-1;k++)
                    for (int ii=0; ii < part[is].npcelx; ii++)
                        for (int jj=0; jj < part[is].npcely; jj++)
                            for (int kk=0; kk < part[is].npcely; kk++){
                                // initialize each particle position and charge. Particle uniform in space
                                part[is].x[counter] = (ii + .5)*(grd->dx/part[is].npcelx) + grd->XN[i][j][k];
                                part[is].y[counter] = (jj + .5)*(grd->dy/part[is].npcely) + grd->YN[i][j][k];
                                part[is].z[counter] = (kk + .5)*(grd->dz/part[is].npcelz) + grd->ZN[i][j][k];
                                // q = charge * statistical weight
                                part[is].q[counter] =  (part[is].qom/fabs(part[is].qom))*(ids[is].rhoc[i][j][k]/part[is].npcel)*(1.0/grd->invVOL);
                                
                                //////////////// Maxwellian ////////////////
                                // u
                                harvest =   rand()/(double)RAND_MAX;
                                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                                harvest =   rand()/(double)RAND_MAX;
                                theta = 2.0*M_PI*harvest;
                                part[is].u[counter] = part[is].u0 + part[is].uth*prob*cos(theta);
                                // check u
                                if (part[is].u[counter] > param->c){
                                    std::cout << "ERROR - u VELOCITY > c !" << std::endl;
                                    exit(EXIT_FAILURE);
                                }
                                // u and w
                                harvest =   rand()/(double)RAND_MAX;
                                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                                harvest =   rand()/(double)RAND_MAX;
                                theta = 2.0*M_PI*harvest;
                                part[is].w[counter] = part[is].w0 + part[is].wth*prob*cos(theta);
                                part[is].v[counter] = part[is].v0 + part[is].vth*prob*sin(theta);
                                // check v and w
                                if (part[is].v[counter] > param->c){
                                    std::cout << "ERROR - v VELOCITY > c !" << std::endl;
                                    exit(EXIT_FAILURE);
                                }
                                if (part[is].w[counter] > param->c){
                                    std::cout << "ERROR - w VELOCITY > c !" << std::endl;
                                    exit(EXIT_FAILURE);
                                }
                                //  update particle counter
                                counter++ ;
                            } // end of one particles initialization
    
    
    } // end of species initialization
    
}

/** initialize uniform electric and magnetic field */
void initUniform(struct parameters* param, struct grid* grd, struct EMfield* field, struct EMfield_aux* field_aux, struct particles* part, struct interpDensSpecies* ids){
    
    // perturbation localized in X
    double pertX = 0.4;
    double xpert, ypert, exp_pert;
    
    // print settings
    std::cout << "*************************************************" << std::endl;
    std::cout << "**  Initialize UNIFORM Plasma Distribution     **" << std::endl;
    std::cout << "*************************************************" << std::endl;
    std::cout << "** B0x = " << param->B0x << std::endl;
    std::cout << "** B0y = " << param->B0y << std::endl;
    std::cout << "** B0z = " << param->B0z << std::endl;
    for (int is=0; is < param->ns; is++){
        std::cout << "** rho species " << is << " = " << param->rhoINIT[is];
    }
    std::cout << "*************************************************" << std::endl;
    
    /////////////////////////////////////////////////
    //////   FIELD AND DENSITY
    /////////////////////////////////////////////////
    // Set the electric field, magnetic field + rhos
    for (int i=0; i < grd->nxn; i++)
        for (int j=0; j < grd->nyn; j++)
            for (int k=0; k < grd->nzn; k++){
                // initialize the density for species
                for (int is=0; is < param->ns; is++){
                        ids[is].rhon[i][j][k] = (FPinterp) param->rhoINIT[is]/param->fourpi;
                }
                //std::cout << "OK" << std::endl;
                // electric field
                // electric field
                field->Ex[i][j][k] =  0.0;
                field->Ey[i][j][k] =  0.0;
                field->Ez[i][j][k] =  0.0;
                field_aux->Exth[i][j][k] =  0.0;
                field_aux->Eyth[i][j][k] =  0.0;
                field_aux->Ezth[i][j][k] =  0.0;
                // Magnetic field
                field->Bxn[i][j][k] = param->B0x;
                field->Byn[i][j][k] = param->B0y;
                field->Bzn[i][j][k] = param->B0z;
                
            }
    // calculate B and rho at centers cells: first argument is on center cell
    interpN2Cfield(field_aux->Bxc, field_aux->Byc, field_aux->Bzc, field->Bxn, field->Byn, field->Bzn, grd);
    

    // interpolate densities species from node
    for (int is=0; is < param->ns; is++){
        interpN2Crho(&ids[is], grd);
    }
    
    /////////////////////////////////////////////////
    //////   PARTICLES
    /////////////////////////////////////////////////
    
    double harvest;
    double prob, theta, sign;
    long long counter;
    
    // loop over the species
    for (int is=0; is < param->ns; is++){
        // set particle counter to zero
        counter=0;
        // set the seed for random number generator equal to species id
        srand(is);
        for (int i=1; i< grd->nxc-1;i++)
            for (int j=1; j< grd->nyc-1;j++)
                for (int k=1; k< grd->nzc-1;k++)
                    for (int ii=0; ii < part[is].npcelx; ii++)
                        for (int jj=0; jj < part[is].npcely; jj++)
                            for (int kk=0; kk < part[is].npcely; kk++){
                                // initialize each particle position and charge. Particle uniform in space
                                part[is].x[counter] = (ii + .5)*(grd->dx/part[is].npcelx) + grd->XN[i][j][k];
                                part[is].y[counter] = (jj + .5)*(grd->dy/part[is].npcely) + grd->YN[i][j][k];
                                part[is].z[counter] = (kk + .5)*(grd->dz/part[is].npcelz) + grd->ZN[i][j][k];
                                // q = charge * statistical weight
                                part[is].q[counter] =  (part[is].qom/fabs(part[is].qom))*(ids[is].rhoc[i][j][k]/part[is].npcel)*(1.0/grd->invVOL);
                                
                                //////////////// Maxwellian ////////////////
                                // u
                                harvest =   rand()/(double)RAND_MAX;
                                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                                harvest =   rand()/(double)RAND_MAX;
                                theta = 2.0*M_PI*harvest;
                                part[is].u[counter] = part[is].u0 + part[is].uth*prob*cos(theta);
                                // check u
                                if (part[is].u[counter] > param->c){
                                    std::cout << "ERROR - u VELOCITY > c !" << std::endl;
                                    exit(EXIT_FAILURE);
                                }
                                // u and w
                                harvest =   rand()/(double)RAND_MAX;
                                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                                harvest =   rand()/(double)RAND_MAX;
                                theta = 2.0*M_PI*harvest;
                                part[is].w[counter] = part[is].w0 + part[is].wth*prob*cos(theta);
                                part[is].v[counter] = part[is].v0 + part[is].vth*prob*sin(theta);
                                // check v and w
                                if (part[is].v[counter] > param->c){
                                    std::cout << "ERROR - v VELOCITY > c !" << std::endl;
                                    exit(EXIT_FAILURE);
                                }
                                if (part[is].w[counter] > param->c){
                                    std::cout << "ERROR - w VELOCITY > c !" << std::endl;
                                    exit(EXIT_FAILURE);
                                }
                                //  update particle counter
                                counter++ ;
                            } // end of one particles initialization
    
    
    } // end of species initialization
    
}
