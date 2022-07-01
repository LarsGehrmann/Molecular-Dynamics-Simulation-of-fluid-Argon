#include "header.h"


double calc_all(size_t particles, size_t dimension, int maxT, double F[particles][3], double Potential[maxT], double r_current[particles][dimension], double v_current[particles][dimension], double boxlength, int number_Of_Bins, double g[number_Of_Bins], int step_Number)
{
	double dr = boxlength / number_Of_Bins;
	double Vol =  boxlength* boxlength* boxlength;
	double sum;
	double Delta_r[3];
	double prefactor;
	double r;
	double r2;
	double r6;
	double h[number_Of_Bins];
	double rs[number_Of_Bins];
	int idx;
	Potential[step_Number] = 0;
	for (int k = 0; k < particles; k++) {		// set Force to zero
		for (int jj = 0; jj < 3; jj++) {
			F[k][jj] = 0;
		}
	}
	for (int j = 0; j < number_Of_Bins; j++) {
		rs[j] = (j+0.5) * dr;
		h[j] = 0;
	}
	for (int j = 0; j < particles; j++) {
		for (int kk = j + 1; kk < particles; kk++) {
			sum = 0;
			for (int iii = 0; iii < 3; iii++) {
			Delta_r[iii] = r_current[j][iii] - r_current[kk][iii];
			Delta_r[iii] -= boxlength * round(Delta_r[iii] / boxlength);
			sum += Delta_r[iii] * Delta_r[iii];
			}
			r2 = 1/sum;
			r6 = r2 * r2 * r2;
			prefactor = 48* r6 * r2* ( r6 - 0.5);	
			for (int ooo = 0; ooo < 3; ooo++) {
				F[j][ooo] += prefactor * Delta_r[ooo];
				F[kk][ooo] -= prefactor * Delta_r[ooo];
			}
			r = sqrt(sum);
			idx = floor(r / dr);
			h[idx] += 1;
			Potential[step_Number] += 4 * r6 * (r6 - 1);
		}
	}
	g[0] = 0;
	for (int i = 1; i < number_Of_Bins; i++) {
		h[i] = h[i] * 2 / particles;
		g[i] = Vol * h[i] / (4 * particles * M_PI * dr * rs[i] * rs[i]);
	}
} // end function
