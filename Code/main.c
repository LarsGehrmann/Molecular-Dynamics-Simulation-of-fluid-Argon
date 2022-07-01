#include "header.h"

void main()
{
	clock_t begin = clock();
	double temperature = 94.4 /119.8 ; // in reduced units
	size_t particles = 125;
	double rho = 1.374;	// in g/cm^3
	double third = 1.0 / 3.0;
	double boxlength = pow(particles * 6.635 * 10 / ( pow(3.405,3) *  rho),third);
	int number_Of_Bins = 100;  // number of bins for correlation function
	printf("boxlength %lf\n", boxlength);
	int maxT = 2000;		// maximum number of timesteps
	delete_all();			// delete previous .txt files
	size_t dimension = 3;
	double dt = 0.001;		// timesteps
	double r_current[particles][3];
	double v_current[particles][3];
	initial(particles, dimension, maxT, r_current, v_current, boxlength,temperature);		// initialize system with periodic lattice and random velocities
	calc_Timesteps(particles, dimension, maxT, r_current, v_current, dt, boxlength,temperature, number_Of_Bins);	// calculate timesteps
	clock_t end = clock();								
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Total time needed: %lfs\n", time_spent);
} // end main