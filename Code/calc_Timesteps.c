#include "header.h"
#include "calc_all.c"
#include "Kinetic_energy.c"
#include "VVs.c"
void calc_Timesteps (size_t particles, size_t dimension, int T_max, double r_current[particles][3], double v_current[particles][3], double dt, double boxlength, double Tdes, int number_Of_Bins)
{	
	FILE *H = NULL;
	FILE *Kin = NULL;
	FILE *Pot = NULL;
	FILE *gs = NULL; 
	FILE *r0x = NULL;
	FILE *r1x = NULL;
	FILE *v0x = NULL;
	FILE *v1x = NULL;
	H = fopen("./data/H.txt", "a");
	Kin = fopen("./data/Kin.txt", "a");
	Pot = fopen("./data/Pot.txt", "a");
	gs = fopen("./data/gs.txt", "a");
	r0x = fopen("./data/r0x.txt", "a");
	r1x = fopen("./data/r1x.txt", "a");
	v0x = fopen("./data/v0x.txt", "a");
	v1x = fopen("./data/v1x.txt", "a");
	double Tcalc;
	double Tratio12;
	double half = 1.0 / 2.0;
	double Hamilton[T_max];
	double Kinetic[T_max];
	double V[T_max];
	double g[number_Of_Bins];
	double F[particles][3];	
	for (int k = 0; k < T_max; k++) {				// for every time step

		Kinetic[k] = Kinetic_energy(particles, dimension, v_current);
		calc_all(particles, dimension, T_max, F, V, r_current, v_current, boxlength, number_Of_Bins, g, k);
		Hamilton[k] = Kinetic[k] + V[k];
	
		fprintf(H, "%lf\n", Hamilton[k]);
		fclose;
		fprintf(Pot, "%lf\n", V[k]);
		fclose;
		fprintf(Kin, "%lf\n", Kinetic[k]);
		fclose;
		/*
		fprintf(r0x, "%lf\n", r_current[0][0]);
		fclose;
		fprintf(r1x, "%lf\n", r_current[1][0]);
		fclose;
		fprintf(v0x, "%lf\n", v_current[0][0]);
		fclose;
		fprintf(v1x, "%lf\n", v_current[1][0]);
		fclose;
		*/
		//printf("step: %d H: %lf Kin: %lf Pot: %lf \n", k, Hamilton[k], Kinetic[k], V[k]);
		
		for (int i = 0; i < number_Of_Bins; i++) {
			fprintf(gs, "%lf ", g[i]);
			fclose;
			}
		fprintf(gs, "\n");
		fclose;
		
		VVs(particles, dimension, r_current, v_current, dt, boxlength, F);
		if (k % 100 == 0) {
			printf("%d\n", k);
		}
		if (k < 1000 && k % 50 == 0) {
			Kinetic[k] = Kinetic_energy(particles, dimension, v_current);
			Tcalc = 2 * Kinetic[k] / (3 * particles);
			Tratio12 = pow((Tdes / Tcalc), half);
			for (int k = 0; k < particles; k++) {
				for (int jj = 0; jj < 3; jj++) {
					v_current[k][jj] = v_current[k][jj] * Tratio12;
				}
			}
		}
	}
	// after time integration
	double T_avg = 0;
	int eqstart = 1500;
	int diff = T_max - eqstart;
	double kin[diff];
	double kin_avg = 0;
	double sigma2 = 0;
	for (int k = eqstart; k < T_max; k++) {
		T_avg += 2 * Kinetic[k] / (3 * particles);
		kin[k-eqstart] = Kinetic[k] / particles;
	}
	T_avg = T_avg / diff;
	for (int j = 0; j < diff; j++) {
		kin_avg += kin[j];
	}
	kin_avg = kin_avg / diff;
	for (int j = 0; j < diff; j++) {
		sigma2 += (kin[j] - kin_avg) * (kin[j] - kin_avg);
	}
	sigma2 = sigma2 / diff;
	printf("Final sigma2: %lf\n", sigma2);
	printf("kin_avg: %lf    T_avg: %lf    sigma2: %lf \n", kin_avg, T_avg, sigma2);
	double cv = 1 / (0.66666 - 4 * particles * sigma2 / (9 * T_avg * T_avg));
	printf("C_V: %lf\n", cv);
}
