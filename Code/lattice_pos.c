void lattice_pos(size_t particles, size_t dimension, double r_current[particles][dimension], double boxlength, int maxT)
{
	double third = 1.0/3.0;
	int after;
	double spaces = pow(particles,third);	// spaces in ONE direction
	spaces = ceil(spaces);
	double a = boxlength/(spaces); // spaces +1 ??
	int counter = 0;
	for (int n = 0; n < spaces; n++) {
		for (int mm = 0; mm < spaces; mm++) {
			for (int ooo = 0; ooo < spaces; ooo++) {
				r_current[counter][0] = n*a;
				r_current[counter][1] = mm*a;
				r_current[counter][2] = ooo*a;
			//	printf("Particle Number %d x:%lf y:%lf z:%lf\n", counter, r_current[counter][0], r_current[counter][1], r_current[counter][2]);
				counter++;
				if (counter == particles)
					goto jump;
				}
		}
	}
jump:
	after = 1;


}
