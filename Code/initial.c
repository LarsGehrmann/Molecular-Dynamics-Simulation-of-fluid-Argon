#include "lattice_pos.c"
#include "ranf.c"

void initial(size_t particles, size_t dimension, int maxT, double r_current[particles][3], double v_current[particles][3], double boxlength, double temperature)
{
	double sumv[3] = {0,0,0};
	double sumv2 = 0;
	double rando;
	double scale;
	lattice_pos(particles, dimension, r_current, boxlength, maxT);
	for (int k = 0; k < particles; k++) {
		for (int jj = 0; jj < 3; jj++) {
			rando = ranf();
			v_current[k][jj] = (rando - 0.5);
			sumv[jj] = sumv[jj] + v_current[k][jj];
			sumv2 = sumv2 + v_current[k][jj]*v_current[k][jj];
		}
	}
	for (int k = 0; k < 3; k++) {
	sumv[k] = sumv[k]/particles;
	}
	sumv2 = sumv2/particles;
	scale = sqrt(3*temperature/sumv2);
	for (int k = 0; k < particles; k++) {
		for (int jj = 0; jj < 3; jj++) {
			v_current[k][jj] = (v_current[k][jj]-sumv[jj])*scale;
		}
	}
}
