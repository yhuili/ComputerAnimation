#include "particle.h"
#include <iostream>

Particle::Particle()
{
	ParticleInit(0.0);
};

Particle::Particle(double mass)
{
	ParticleInit(mass);
};

void Particle::ParticleInit(double mass)
{
	p_mass = mass;

	p_position.x = 0.0;
	p_position.y = 0.0;
	p_position.z = 0.0;

	p_velocity.x = 0.0;
	p_velocity.y = 0.0;
	p_velocity.z = 0.0;

	p_acceleration.x = 0.0;
	p_acceleration.y = 0.0;
	p_acceleration.z = 0.0;
};

Particle::~Particle(){};

ParticleSystem::ParticleSystem(int num)
{
	particle_num = num;
	link_length = 1.0 / (particle_num - 1);
	double mass = 1.0 / particle_num;
	Particle p(mass);
	particle_list.resize(particle_num, p);
	particle_list[0].p_mass = 0.0; // set the first particle mass to 0
	// initialize all variables
	integrator = 'E'; // 'R' or 'E'
	gravity.x = 0.0;
	gravity.y = -1.0;
	gravity.z = 0.0;
	k_damp = 0.5;
	dt = 0.001;
	for (int i = 0; i < 4; i++)
	{
		link_dc_dq[i] = 0.0;
		link_dcDot_dq[i] = 0.0;
	}
	for (int i = 0; i < 2; i++)
	{
		ring_dc_dq[i] = 0.0;
		ring_dcDot_dq[i] = 0.0;
	}
	row = particle_num; // particle_num = 11
	col = (particle_num - 1) * 2;
	alpha = 5.0; // 2b
	beta = (alpha * alpha) / 4.0; // b*b
	// ring variables
	inner_radius = 0.0075;
	outter_radius = 0.5;
	num_minor = 30;
	num_major = 30;
	ring_position.x = 0.0;
	ring_position.y = -0.5;
	ring_position.z = 0.0;
	error_value.open("error.txt");
	frame = 0;
	up_external_force = 0.0;
	down_external_force = 0.0;
	left_external_force = 0.0;
	right_external_force = 0.0;

	// set the position of all particles to 3 o'clock
	for (int i = 0; i < (particle_num - 1) / 2; i++)
	{
		particle_list[i].p_position.x = 0.0;
		particle_list[i].p_position.y = particle_list[i].p_position.y - link_length * i;
		particle_list[i].p_position.z = 0.0;
	}

	for (int i = (particle_num - 1) / 2; i < particle_num; i++)
	{
		particle_list[i].p_position.x = particle_list[i].p_position.x + link_length * (i - (particle_num - 1) / 2);
		particle_list[i].p_position.y = -0.5;
		particle_list[i].p_position.z = 0.0;
	}
	// allocate memory to all matrix
	mass_matrix = new double[col * col];
	memset(mass_matrix, 0, sizeof(double) * col * col);
	
	gradientC_matrix = new double[row * col];
	memset(gradientC_matrix, 0, sizeof(double) * row * col);
	
	gradientCDot_matrix = new double[row * col];
	memset(gradientCDot_matrix, 0, sizeof(double) * row * col);
	
	transposed_gradientC_matrix = new double[col * row];
	memset(transposed_gradientC_matrix, 0, sizeof(double) * col * row);
	
	int new_dimension = col + row;
	left_matrix = new double[new_dimension * new_dimension];
	memset(left_matrix, 0, sizeof(double) * new_dimension * new_dimension);
	
	external_force_matrix = new double[col];
	qDot_matrix = new double[col];
	memset(qDot_matrix, 0, sizeof(double) * col);
	
	gradientCqDotDot_matrix_1 = new double[row];
	memset(gradientCqDotDot_matrix_1, 0, sizeof(double) * row);
	
	gradientCqDotDot_matrix_2 = new double[row];
	memset(gradientCqDotDot_matrix_2, 0, sizeof(double) * row);
	
	gradientCqDotDot_matrix_3 = new double[row];
	memset(gradientCqDotDot_matrix_3, 0, sizeof(double) * row);
	
	right_matrix = new double[new_dimension];
	memset(right_matrix, 0, sizeof(double) * new_dimension);
	
	x = gsl_vector_alloc(new_dimension);
	gsl_vector_set_zero(x);
};

ParticleSystem::~ParticleSystem()
{
	delete [] mass_matrix;
	delete [] gradientC_matrix;
	delete [] gradientCDot_matrix;
	delete [] transposed_gradientC_matrix;
	delete [] left_matrix;
	delete [] external_force_matrix;
	delete [] qDot_matrix;
	delete [] gradientCqDotDot_matrix_1;
	delete [] gradientCqDotDot_matrix_2;
	delete [] gradientCqDotDot_matrix_3;
	delete [] right_matrix;
	gsl_vector_free(x);
};

double ParticleSystem::LinkConstraint(int index)
{
	// start from the third particle in the particle_list (index = 2), particle_list[0].p_position = (0, 0, 0)T
	double x = particle_list[index - 1].p_position.x - particle_list[index].p_position.x;
	double y = particle_list[index - 1].p_position.y - particle_list[index].p_position.y;
	double c  = x * x + y * y - link_length * link_length;
	return c;
};

double ParticleSystem::RingConstraint(int index)
{
	// index = last particle# (n = 10), c = ||xn - p|| - 0.5, p = (0, -0.5, 0)T
	double c = particle_list[index].p_position.x * particle_list[index].p_position.x 
		+ particle_list[index].p_position.y * particle_list[index].p_position.y + particle_list[index].p_position.y;
	return c;
};

void ParticleSystem::LinkGradient(int index)
{
	// start from the third particle in the particle_list (index = 2) to the end (index = 10)
	link_dc_dq[0] = 2 * (particle_list[index - 1].p_position.x - particle_list[index].p_position.x);
	link_dc_dq[1] = 2 * (particle_list[index - 1].p_position.y - particle_list[index].p_position.y);
	link_dc_dq[2] = 2 * (particle_list[index].p_position.x - particle_list[index - 1].p_position.x);
	link_dc_dq[3] = 2 * (particle_list[index].p_position.y - particle_list[index - 1].p_position.y);
};

void ParticleSystem::RingGradient(int index)
{
	// index = last particle# (n = 10), c' = (2xn, 2yn + 1)
	ring_dc_dq[0] = 2 * particle_list[index].p_position.x;
	ring_dc_dq[1] = 2 * particle_list[index].p_position.y + 1;
};

void ParticleSystem::ComputeLinkCDot(int index)
{
	// start from the third particle in the particle_list (index = 2) to the end (index = 10)
	link_dcDot_dq[0] = 2 * (particle_list[index - 1].p_velocity.x - particle_list[index].p_velocity.x);
	link_dcDot_dq[1] = 2 * (particle_list[index - 1].p_velocity.y - particle_list[index].p_velocity.y);
	link_dcDot_dq[2] = 2 * (particle_list[index].p_velocity.x - particle_list[index - 1].p_velocity.x);
	link_dcDot_dq[3] = 2 * (particle_list[index].p_velocity.y - particle_list[index - 1].p_velocity.y);
};

void ParticleSystem::ComputeRingCDot(int index)
{
	// index = last particle# (n = 10)
	ring_dcDot_dq[0] = 2 * particle_list[index].p_velocity.x;
	ring_dcDot_dq[1] = 2 * particle_list[index].p_velocity.y;
};

void ParticleSystem::MassMatrix()
{
	double m = particle_list[1].p_mass;
	int i, j;
	// mass matrix is a col * col matrix
	for (i = 0, j = 0; i < col; i++, j++)
	{
		mass_matrix[i * col + j] = m;
	}
};

void ParticleSystem::GradientCMatrix()
{
	int index = 0;
	int index_n = particle_num - 1; // the last particle index in the particle list vector
	int size = row * col; // gradientC matrix is a row * col matrix
	// the second particle
	gradientC_matrix[0] = 2 * particle_list[1].p_position.x;
	gradientC_matrix[1] = 2 * particle_list[1].p_position.y;
	// from the third to the last
	for (int i = 2; i < particle_num; i++)
	{
		index = (i - 1) * col + 2 * (i - 2);
		LinkGradient(i);
		gradientC_matrix[index] = link_dc_dq[0];
		gradientC_matrix[index + 1] = link_dc_dq[1];
		gradientC_matrix[index + 2] = link_dc_dq[2];
		gradientC_matrix[index + 3] = link_dc_dq[3];
	}
	// ring constraint to the last particle
	RingGradient(index_n);
	gradientC_matrix[size - 2] = ring_dc_dq[0];
	gradientC_matrix[size - 1] = ring_dc_dq[1];
};

void ParticleSystem::GradientCDotMatrix()
{
	int index = 0;
	int index_n = particle_num - 1; // the last particle index in the particle list vector
	int size = row * col; // gradientC matrix is a row * col matrix
	// the first particle
	gradientCDot_matrix[0] = 2 * particle_list[1].p_velocity.x;
	gradientCDot_matrix[1] = 2 * particle_list[1].p_velocity.y;
	// from the second to the last
	for (int i = 2; i < particle_num; i++)
	{
		index = (i - 1) * col + 2 * (i - 2);
		ComputeLinkCDot(i);
		gradientCDot_matrix[index] = link_dcDot_dq[0];
		gradientCDot_matrix[index + 1] = link_dcDot_dq[1];
		gradientCDot_matrix[index + 2] = link_dcDot_dq[2];
		gradientCDot_matrix[index + 3] = link_dcDot_dq[3];
	}
	// ring constraint to the last particle
	ComputeRingCDot(index_n);
	gradientCDot_matrix[size - 2] = ring_dcDot_dq[0];
	gradientCDot_matrix[size - 1] = ring_dcDot_dq[1];
};

void ParticleSystem::TransposedGradientCMatrix()
{
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			transposed_gradientC_matrix[j * row + i] = gradientC_matrix[i * col + j];
		}
	}
};

void ParticleSystem::LeftMatrix()
{
	int i, j, k;
	int new_dimension = row + col;
	// insert mass matrix to left matrix
	MassMatrix();
	for (i = 0; i < col; i++)
	{
		for (j = 0; j < col; j++)
		{
			left_matrix[i * new_dimension + j] = mass_matrix[i * col + j];
		}
	}
	// insert gradientC matrix to left matrix
	GradientCMatrix();
	for (i = col, k = 0; i < new_dimension, k < row; i++, k++)
	{
		for (j = 0; j < col; j++)
		{
			left_matrix[i * new_dimension + j] = gradientC_matrix[k * col + j];
		}
	}
	// insert transposed gradientC matrix to left matrix
	TransposedGradientCMatrix();
	for (i = 0; i < col; i++)
	{
		for (j = col, k = 0; j < new_dimension, k < row; j++, k++)
		{
			left_matrix[i * new_dimension + j] = transposed_gradientC_matrix[i * row + k];
		}
	}
};

void ParticleSystem::ExternalForceMatrix()
{
	memset(external_force_matrix, 0, sizeof(double) * col);
	// gravity force on element y
	double m = particle_list[1].p_mass;
	int i, j;
	for (i = 1; i < col; i += 2)
	{
		external_force_matrix[i] = external_force_matrix[i] + gravity.y * m; //external_force_matrix[i] = gravity.y * m;
	}

	//// damping force on x and y
	//for (i = 0, j = 1; i < col, j < particle_num; i += 2, j++)
	//{
	//	external_force_matrix[i] = external_force_matrix[i] + particle_list[j].p_velocity.x * (-1.0) * k_damp;
	//	external_force_matrix[i + 1] = external_force_matrix[i + 1] + particle_list[j].p_velocity.y * (-1.0) *k_damp;
	//}

	// UI control force
	for (i = 0; i < col; i += 2)
	{
		external_force_matrix[i] = external_force_matrix[i] - left_external_force * m;
		external_force_matrix[i] = external_force_matrix[i] + right_external_force * m;
		external_force_matrix[i + 1] = external_force_matrix[i + 1] + up_external_force * m;
		external_force_matrix[i + 1] = external_force_matrix[i + 1] - down_external_force * m;
	}
};

void ParticleSystem::MatrixMultiplication(double *matrix_top, double *matrix_new, double *matrix_out)
{
	// allocate memory
	gsl_matrix *gsl_top = gsl_matrix_alloc(row, col);
	gsl_vector *gsl_new = gsl_vector_alloc(col);
	gsl_vector *gsl_out = gsl_vector_alloc(row);

	// set matrix top
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			gsl_matrix_set(gsl_top, i, j, matrix_top[i * col + j]);
		}
	}
	// set matrix new
	for (int i = 0; i < col; i++)
	{
		gsl_vector_set(gsl_new, i, matrix_new[i]);
	}
	// get matrix out by multiplication
	gsl_blas_dgemv(CblasNoTrans, 1, gsl_top, gsl_new, 0.0, gsl_out);
	for (int i = 0; i < row; i++)
	{
		matrix_out[i] = gsl_vector_get(gsl_out, i);
	}

	// free memory
	gsl_matrix_free(gsl_top);
	gsl_vector_free(gsl_new);
	gsl_vector_free(gsl_out);
};

void ParticleSystem::RightMatrix()
{
	int i, j;
	int new_dimension = row + col;
	int index_n = particle_num - 1; // the last particle index in the particle list vector
	
	ExternalForceMatrix();
	// insert external force to right matrix
	for (i = 0; i < col; i++)
	{
		right_matrix[i] = external_force_matrix[i];
	}
	// build matrix q'
	for (i = 0, j = 1; i < col, j < particle_num; i += 2, j++)
	{
		qDot_matrix[i] = particle_list[j].p_velocity.x;
		qDot_matrix[i + 1] = particle_list[j].p_velocity.y;
	}
	// build matrix dC' / dq
	GradientCDotMatrix();
	// compute -(dC' / dq)q'
	MatrixMultiplication(gradientCDot_matrix, qDot_matrix, gradientCqDotDot_matrix_1);
	for (i = 0; i < row; i++)
	{
		gradientCqDotDot_matrix_1[i] = (-1.0) * gradientCqDotDot_matrix_1[i]; 
	}
	// compute (dC / dq)q'
	MatrixMultiplication(gradientC_matrix, qDot_matrix, gradientCqDotDot_matrix_2);
	// compute C
	gradientCqDotDot_matrix_3[0] = particle_list[1].p_position.x * particle_list[1].p_position.x + particle_list[1].p_position.y * particle_list[1].p_position.y - link_length * link_length;
	for (i = 2; i < particle_num; i++)
	{
		gradientCqDotDot_matrix_3[i - 1] = LinkConstraint(i);
	}
	gradientCqDotDot_matrix_3[row - 1] = RingConstraint(index_n);
	// insert -(dC' / dq)q' - alpha * (dC / dq)q' - beta * C to right matrix
	for (i = 0; i < row; i++)
	{
		right_matrix[i + col] = gradientCqDotDot_matrix_1[i] - alpha * gradientCqDotDot_matrix_2[i] - beta * gradientCqDotDot_matrix_3[i];
	}
};

void ParticleSystem::ComputeAcceleration()
{
	int i, j;
	int new_dimension = row + col;
	LeftMatrix(); // get left matrix, (col + row) * (col + row)
	RightMatrix(); // get right matrix, (col + row) * 1
	// allocate memory for V, S, U
	gsl_matrix *gsl_v = gsl_matrix_alloc(new_dimension, new_dimension);
	gsl_vector *gsl_s = gsl_vector_alloc(new_dimension);
	gsl_vector *gsl_work = gsl_vector_alloc(new_dimension);

	gsl_matrix_view m = gsl_matrix_view_array(left_matrix, new_dimension, new_dimension);
	gsl_vector_view b = gsl_vector_view_array(right_matrix, new_dimension);

	// A = U S transposedV
	gsl_linalg_SV_decomp(&m.matrix, gsl_v, gsl_s, gsl_work);
	
	// filter the singular values, truncate to 0 any singular value that is less than eps * largest singular value. 
	// The largest element is S0, and eps = 1E-6
	for (i = 0; i < new_dimension; i++)
	{
		if (fabs(gsl_vector_get(gsl_s, i) / gsl_vector_get(gsl_s, 0)) < 1E-6)
		{
			gsl_vector_set(gsl_s, i, 0.0);
		}
	}
	// solve the equation
	gsl_linalg_SV_solve(&m.matrix, gsl_v, gsl_s, &b.vector, x);

	// only need the first col elements from x vector
	for (i = 0, j = 1; i < col, j < particle_num; i += 2, j++)
	{
		particle_list[j].p_acceleration.x = gsl_vector_get(x, i);
		particle_list[j].p_acceleration.y = gsl_vector_get(x, i + 1);
	}

	// free memory
	gsl_matrix_free(gsl_v);
	gsl_vector_free(gsl_s);
	gsl_vector_free(gsl_work);

	// reset result
	memset(left_matrix, 0, sizeof(double) * new_dimension * new_dimension);
	memset(right_matrix, 0, sizeof(double) * new_dimension);
	gsl_vector_set_zero(x);
};

void ParticleSystem::PrintError()
{
	double error = 0.0;
	for (int i = 0; i < particle_num; ++i)
	{
		// constraint may not always be zero
		error += gradientCqDotDot_matrix_3[i];
		}
	error_value << error << endl;
	frame++;
	if (frame >= 300)
	{
		error_value.close();
	}
}

void ParticleSystem::glDrawCylinder(float c_radius, int c_subdivisions)
{
	float vx;
	float vy;
	float vz;
	float v;
	float ax;

	GLUquadricObj *quadratic = gluNewQuadric();
	gluQuadricNormals(quadratic, GLU_SMOOTH);
	for (int i = 0, j = 1; j < particle_num; i++, j++)
	{
		vx = particle_list[j].p_position.x - particle_list[i].p_position.x;
		vy = particle_list[j].p_position.y - particle_list[i].p_position.y;
		vz = particle_list[j].p_position.z - particle_list[i].p_position.z;
		v = sqrt(vx * vx + vy * vy + vz * vz);

		if (fabs(vz) < 1.0e-3) 
		{
			ax = 57.2957795 * acos(vx / v); // rotation angle in x-y plane
			if (vy <= 0.0)
				ax = -ax;
		}
		else 
		{
			ax = 57.2957795 * acos(vz / v); // rotation angle
			if (vz <= 0.0)
				ax = -ax;
		}

		float rx = -vy * vz;
		float ry = vx * vz;
		glPushMatrix();

		// Draw the cylinder body
		glTranslated(particle_list[i].p_position.x, particle_list[i].p_position.y, particle_list[i].p_position.z);
		if (fabs(vz) < 1.0e-3) 
		{
			glRotated(90.0, 0, 1, 0.0); // Rotate & align with x axis
			glRotated(ax, -1.0, 0.0, 0.0); // Rotate to point 2 in x-y plane
		}
		else 
		{
			glRotated(ax, rx, ry, 0.0); // Rotate about rotation vector
		}
		gluQuadricOrientation(quadratic, GLU_OUTSIDE);
		gluCylinder(quadratic, c_radius, c_radius, v, c_subdivisions, 1);

		// Draw the first cap
		gluQuadricOrientation(quadratic, GLU_INSIDE);
		gluDisk(quadratic, 0.0, c_radius, c_subdivisions, 1);
		glTranslated(0, 0, v);

		// Draw the second cap
		gluQuadricOrientation(quadratic, GLU_OUTSIDE);
		gluDisk(quadratic, 0.0, c_radius, c_subdivisions, 1);
		glPopMatrix();
	}
	gluDeleteQuadric(quadratic);
};

void ParticleSystem::glDrawSphere(float s_radius, int s_subdivisions)
{
	for (int i = 0; i < particle_num; i++)
	{
		GLUquadricObj *quadratic = gluNewQuadric();
		gluQuadricNormals(quadratic, GLU_SMOOTH);
		glPushMatrix(); // (NEW) create new matrix
		glTranslated(particle_list[i].p_position.x, particle_list[i].p_position.y, particle_list[i].p_position.z);
		gluSphere(quadratic, s_radius, s_subdivisions, s_subdivisions);
		glPopMatrix();
		gluDeleteQuadric(quadratic);
	}
	
}

void ParticleSystem::glDrawRing()
{
	point normal;
	double major_step = 2.0f * pi / num_major;
	double minor_step = 2.0f * pi / num_minor;

	for (int i = 0; i < num_major; i++)
	{
		double a0 = i * major_step;
		double a1 = a0 + major_step;
		float x0 = (float) cos(a0);
		float y0 = (float) sin(a0);
		float x1 = (float) cos(a1);
		float y1 = (float) sin(a1);
		glBegin(GL_TRIANGLE_STRIP);
		
		for (int j = 0; j <= num_minor; j++)
		{
			double b = j * minor_step;
			float c = (float) cos(b);
			float r = inner_radius * c + outter_radius;
			float z = inner_radius * (float) sin(b);

			glTexCoord2f((float)i / (float)(num_major), (float)(j) / (float)(num_minor));
			normal.x = x0 * c;
			normal.y = y0 * c;
			normal.z = z / inner_radius;

			// normalize normal
			double length = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
			normal.x = normal.x / length;
			normal.y = normal.y / length;
			normal.z = normal.z / length;

			glNormal3f(normal.x, normal.y, normal.z);
			glVertex3f(x0 * r + ring_position.x, y0 * r + ring_position.y, z + ring_position.z);

			glTexCoord2f((float)(i + 1) / (float)(num_major), (float)(j) / (float)(num_minor));
			normal.x = x1 * c;
			normal.y = y1 * c;
			normal.z = z / inner_radius;

			length = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
			normal.x = normal.x / length;
			normal.y = normal.y / length;
			normal.z = normal.z / length;

			glNormal3f(normal.x, normal.y, normal.z);
			glVertex3f(x1 * r + ring_position.x, y1 * r + ring_position.y, z + ring_position.z);
		}
		glEnd();
	}
}

void ParticleSystem::glRender(float s_radius, int s_subdivisions, float c_radius, int c_subdivisions)
{
	glDrawCylinder(c_radius, c_subdivisions);
	glDrawSphere(s_radius,s_subdivisions);
	glDrawRing();
};


