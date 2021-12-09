#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
/* #define SiO2_MASS 60.0843
#define MASS_CONVERTION_FACTOR 9649.0 */
#define PRESSURE 99.8*1e12   // ag/(mm*(ms)^2)
#define KB 1.380649*1e-2 //	[ag*(mm)^2/(ms)^2])// 8.617333262145*1e-5 // eV/K, 	
#define dims 1

/* Simulation parameters */
const double dt = 0.001; /* ms = milli sec = 10^-3 sec */
const double tau = 147.3*1e-3; /* ms = milli sec */
const double friction_coefficient = 1/tau;
const double omega_0 = 2.0*M_PI * 3.1; /* 1/ms */
const double T = 297.0; /* K */
const long int sim_len = (long int)(2000)/((long int)(dt*1e3));
const long int n_timesteps = 3000*sim_len; /* X * [production group steps] */



/* Properties of Brownian particle */
typedef struct Brownian_particle{
    double r[dims];
    double v[dims];
    double a[dims];
    double m;
    double m_inv;
    double d;
    double rho;
}Brownian_particle;



void velocity_verlet_algorithm(double **, double **, double, 
                            Brownian_particle, double, 
                            gsl_rng *, gsl_rng *, gsl_rng *);



int main(){

    const double c_0 = 1.0/exp(friction_coefficient * dt);

    /* Initiate particle */
    struct Brownian_particle p1;
    for(int ix = 0; ix < dims; ++ix){
        p1.r[ix] = 0;
    }
    for(int ix = 0; ix < dims; ++ix){
        p1.v[ix] = 0;
    }
    for(int ix = 0; ix < dims; ++ix){
        p1.a[ix] = 0;
    }
    p1.m =  30134022.3516; /* ag = atto gram = 10^-18g */  /* SiO2_MASS/MASS_CONVERTION_FACTOR; */
    p1.d = 0.5 * 2.79 * 1e-3; // [r] = mm
    p1.rho = 2.65*1e15; // ag/(mm)^3
    p1.m_inv = 1.0/p1.m;

    /* Random number generator */
    const gsl_rng_type *randType1 = gsl_rng_default;
    unsigned long int seed1 = 16384;
    gsl_rng *rand_gen1 = gsl_rng_alloc(randType1);
    gsl_rng_set( rand_gen1, seed1);
    gsl_ran_ugaussian( rand_gen1 );


    const gsl_rng_type *randType2 = gsl_rng_default;
    unsigned long int seed2 = 947223;
    gsl_rng *rand_gen2 = gsl_rng_alloc(randType2);
    gsl_rng_set( rand_gen2, seed2);
    gsl_ran_ugaussian( rand_gen2 );

    const gsl_rng_type *randType3 = gsl_rng_default;
    unsigned long int seed3 = 128;
    gsl_rng *rand_gen3 = gsl_rng_alloc(randType3);
    gsl_rng_set( rand_gen3, seed3);
    gsl_ran_ugaussian( rand_gen3 );

    /* Just for testing the random numbers */
    /* run the program with the command " ./task1 > ugaussian_results.dat "
        to save prinf(...) output to file */
    /*
    long int iter = 10000;
    for( long int ix = 0; ix < iter; ++ix ){
            printf("%lf\n", gsl_ran_ugaussian( rand_gen ));
    }
    */

    /* Malloc to store trajectory in phase space */

    double * v_phase_entries = (double *)malloc( n_timesteps * dims 
    * sizeof(double) );
    double **v_phase = (double **)malloc( n_timesteps * sizeof(double *) );
    for(long int row = 0, jump = 0; row < n_timesteps; ++row, jump += dims)
    {
        v_phase[row] = v_phase_entries + jump;
    }
     for(long int tx = 0; tx < n_timesteps; ++tx)
    {
        for(int dim = 0; dim < dims; ++dim)
            v_phase[tx][dim] = 0;
    }


    double * r_phase_entries = (double *)malloc(
        n_timesteps * dims * sizeof(double) );

    double **r_phase = (double **)malloc(
        n_timesteps * sizeof(double *) );
    for(long int row = 0, jump = 0; row < n_timesteps; ++row, jump += dims)
    {
        r_phase[row] = r_phase_entries + jump;
    }
     for(long int tx = 0; tx < n_timesteps; ++tx)
    {
        for(int dim = 0; dim < dims; ++dim)
            r_phase[tx][dim] = 0;
    }

    /* Relaxation */
    velocity_verlet_algorithm( r_phase, v_phase, c_0, p1, dt, 
    rand_gen1, rand_gen1, rand_gen1);
    velocity_verlet_algorithm( r_phase, v_phase, c_0, p1, dt, 
    rand_gen1, rand_gen1, rand_gen1);

    /* Run velocity verlet on a particle */
    velocity_verlet_algorithm( r_phase, v_phase, c_0, p1, dt, 
    rand_gen1, rand_gen1, rand_gen1);

    printf("%lf,%li", dt, sim_len);
    for(long int tx = 0; tx < n_timesteps; ++tx)
    {
        printf("\n%lf,%lf", 
        r_phase[tx][0], v_phase[tx][0]);
    }

    free(r_phase_entries);
    free(r_phase);
    free(v_phase_entries);
    free(v_phase);
    gsl_rng_free (rand_gen1);
    gsl_rng_free (rand_gen2);
    gsl_rng_free (rand_gen3);

}



void velocity_verlet_algorithm( double ** r_phase, double ** v_phase, 
            double c_0, Brownian_particle p, double dt, 
            gsl_rng *rand_gen1, gsl_rng *rand_gen2, gsl_rng *xi)
{   
    
    double gauss_rand;
    double v_th = sqrt(KB*T*p.m_inv);
    for(long int tx = 0; tx < n_timesteps; ++tx){
        
        /* half-step */
        for(int dim = 0; dim < dims; ++dim)
        {   
            gauss_rand = gsl_ran_ugaussian( rand_gen1 );

            p.v[dim] = 0.5*p.a[dim]*dt 
                        + sqrt(c_0)*p.v[dim] 
                        + v_th*sqrt(1 - c_0)*gauss_rand;

            p.r[dim] = p.r[dim] + p.v[dim]*dt;
        }

        /* Calculate forces */
        for(int dim = 0; dim < dims; ++dim)
        {
            gauss_rand = gsl_ran_ugaussian( rand_gen1 );
            p.a[dim] = -(omega_0*omega_0)*p.r[dim] 
                        - friction_coefficient*p.v[dim]
                        + gauss_rand;
        }

        /* Full-step */
        for(int dim = 0; dim < dims; ++dim)
        {   
            gauss_rand = gsl_ran_ugaussian( rand_gen1 );

            p.v[dim] = 0.5*sqrt(c_0)*p.a[dim]*dt 
                        + sqrt(c_0)*p.v[dim] 
                        + v_th*sqrt(1 - c_0)*gauss_rand;

        }

        /* Save trajectory */
        for(int dim = 0; dim < dims; ++dim)
            r_phase[tx][dim] = p.r[dim];
        
        for(int dim = 0; dim < dims; ++dim)
            v_phase[tx][dim] = p.v[dim];
    }
}