#include "Force_analytical.c"

void  VVs(size_t particles, size_t dimension, double r_current[particles][dimension], double v_current[particles][dimension],double dt, double boxlength, double Force[particles][3])
{
	// (1) calculate v_half
	double v_half[particles][3];
	for (int k = 0; k < particles; k++) {		// for every particle
		for (int kk = 0; kk < 3; kk++) {
			v_half[k][kk] = v_current[k][kk] + 0.5*Force[k][kk]*dt;
		}
	}
	// (2) calculate r_next
	for (int k = 0; k < particles; k++) {		// for every particle
		for (int kk = 0; kk < 3; kk++) {
			r_current[k][kk] += v_half[k][kk]*dt;
		}
	}
	// (3) calculate a_next
	Force_analytical(particles, dimension, Force, r_current, boxlength);

	// (4) calculate v_next
	for (int k = 0; k < particles; k++) {		// for every particle
		for (int kk = 0; kk < 3; kk++) {
			v_current[k][kk] = v_half[k][kk] + 0.5*Force[k][kk]*dt;
		}
	}	
} // end function
