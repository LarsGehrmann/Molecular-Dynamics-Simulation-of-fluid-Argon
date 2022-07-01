#ifndef header
# define header
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <dirent.h>
#include <string.h>
#include <time.h>
#define m 1 // mass in reduced units
#define sigma 1 // in reduced units
#define eps 1 // in reduced units
#include "delete_all.c"
#include "initial.c"
#include "calc_Timesteps.c"

double ranf();
double norm(double r[3]);
void lattice_pos(size_t particles, size_t dimension, double r_current[particles][dimension], double boxlength, int maxT);
void subtract(double Delta_r[3], double r1[3], double r2[3]);
double V_total(size_t particles, size_t dimension, double rmatrix[particles][dimension],double boxlength);
void Force_analytical(size_t particles, size_t dimension, double F[particles][3], double rmatrix[particles][dimension],  double boxlength);
double Kinetic_energy(size_t particles, size_t dimension, double v_current[particles][3]);
void  VVs(size_t particles, size_t dimension, double r_current[particles][dimension], double v_current[particles][dimension], double dt, double boxlength, double Force[particles][3]);
void delete_all();
void calc_Timesteps(size_t particles, size_t dimension, int T_max, double r_current[particles][3], double v_current[particles][3], double dt, double boxlength, double Tdes, int number_Of_Bins);
void initial(size_t particles, size_t dimension, int maxT, double r_current[particles][3], double v_current[particles][3], double boxlength, double temperature);
void calc_RDF(size_t particles, size_t dimension, double r_current[particles][3], double boxlength, int number_Of_Bins, double g[number_Of_Bins]);
double calc_all(size_t particles, size_t dimension, int maxT, double F[particles][3], double Potential[maxT], double r_current[particles][dimension], double v_current[particles][dimension], double boxlength, int number_Of_Bins, double g[number_Of_Bins], int step_Number);





#endif