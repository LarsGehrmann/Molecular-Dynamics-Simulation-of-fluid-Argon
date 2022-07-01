double Kinetic_energy(size_t particles, size_t dimension, double v_current[particles][3]) {
	double Kinetic = 0;
	
	for (int k = 0; k < particles; k++) {
		for (int jj = 0; jj < 3; jj++) {
			Kinetic += 0.5 * v_current[k][jj] * v_current[k][jj];
		}
	}
	return Kinetic;
/*
	double v_abs[particles];
double v_help[3];
	for (int k = 0; k < particles; k++) {
		for (int jj = 0; jj < 3; jj++) {
			v_help[jj] = v_current[k][jj];
		}
		v_abs[k] = norm(v_help);
		Kinetic += 0.5*m*v_abs[k]*v_abs[k];
	}
	return Kinetic;
	*/
} // end function
