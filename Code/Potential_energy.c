#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double Potential_energy(size_t particles, size_t dimension,double r_current[particles][dimension], double boxlength, int Particle_Number)
{
	double Delta_r[3];
	double prefactor;
	double r;
	double r2;
	double r6;
	double sum = 0;
	double Potential = 0;
	double r_current_part[3];
	
	for (int j = 0; j < 3; j++) {
		r_current_part[j] = r_current[Particle_Number][j];
	}
	for (int j = 0; j < particles; j++) {
		if (j != Particle_Number) {
			sum = 0;
			for (int ii = 0; ii < 3; ii++) {
					Delta_r[ii] = r_current_part[ii] - r_current[j][ii];
					Delta_r[ii] -= boxlength * round(Delta_r[ii] / boxlength);
					sum += Delta_r[ii] * Delta_r[ii];
				}
				r2 = 1 / sum;
				r6 = r2 * r2 * r2;
				Potential += 4 * r6 * (r6 - 1);
			}
			else
				;
		}
	return Potential;

} // end function
