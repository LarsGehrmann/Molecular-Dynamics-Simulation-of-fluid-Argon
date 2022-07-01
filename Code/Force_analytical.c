#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void Force_analytical(size_t particles, size_t dimension, double F[particles][3], double r_current[particles][dimension],  double boxlength)
{
	double Delta_r[3];
	double r2;
	double r6;
	double prefactor;
	double sum;
	
	for (int k = 0; k < particles; k++) {		// set Force to zero
		for (int jj = 0; jj < 3; jj++) {
			F[k][jj] = 0;
		}
	}
	
	for (int j = 0; j < particles; j++) {
		for (int kk = j + 1; kk < particles; kk++) {
			sum = 0;
			for (int iii = 0; iii < 3; iii++) {
				Delta_r[iii] = r_current[j][iii] - r_current[kk][iii];
				Delta_r[iii] -= boxlength * round(Delta_r[iii] / boxlength);
				sum += Delta_r[iii] * Delta_r[iii];
			}
			r2 = 1 / sum;
			r6 = r2 * r2 * r2;
			prefactor = 48 * r6 * r2 * (r6 - 0.5);
			for (int ooo = 0; ooo < 3; ooo++) {
				F[j][ooo] += prefactor * Delta_r[ooo];
				F[kk][ooo] -= prefactor * Delta_r[ooo];
			}
		}
	}
} // end function
