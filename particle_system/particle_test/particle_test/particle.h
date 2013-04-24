#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include "jello.h"
#include <vector>
#include <fstream>

using namespace std;

#ifndef _PARTICLE_H_
#define _PARTICLE_H_
class Particle
{
public:
	double p_mass; // the mass of the particle point
	point p_position; // position of the particle point
	point p_velocity; // velocities of the particle point
	point p_acceleration; // acceleration of the particle point

	Particle(); // particle constructor
	Particle(double mass);
	~Particle(); // particle destructor
private:
	void ParticleInit(double mass);
};
#endif

#ifndef _PARTICLESYSTEM_H_
#define _PARTICLESYSTEM_H_
class ParticleSystem
{
public:
	// member variable
	int particle_num; // number of particles, public
	vector<Particle> particle_list; // the list of particles in the system
	char integrator; // Euler or RK4
	double dt; // simulation clock
	ofstream error_value;
	int frame;
	// UI control
	double left_external_force;
	double right_external_force;
	double up_external_force;
	double down_external_force;

	// member function
	// particle system constructor
	//ParticleSystem();
	ParticleSystem(int num);
	~ParticleSystem(); // particle system destructor
	// solve the equation using gsl to get the acceleration
	void ComputeAcceleration();
	// particle system renderer
	void glRender(float s_radius, int s_subdivisions, float c_radius, int c_subdivisions);
	// write the error to txt file
	void PrintError();

private:
	double link_length; // length of the link between each two particles
	point gravity; // gravity force
	double k_damp; // damping constant k
	double link_dc_dq[4];
	double ring_dc_dq[2];
	double link_dcDot_dq[4];
	double ring_dcDot_dq[2];
	int row; // the row # of the matrix
	int col; // the col # of the matrix
	//  build the left matrix
	double *mass_matrix; // M
	double *gradientC_matrix; // dC / dq
	double *gradientCDot_matrix; // dC' / dq
	double *transposed_gradientC_matrix; // (dC / dq)T
	double *left_matrix;
	// build the right matrix
	double *external_force_matrix; // f(t)
	double *qDot_matrix; // q'
	double *gradientCqDotDot_matrix_1; // -(dC' / dq)q'
	double *gradientCqDotDot_matrix_2; // (dC / dq)
	double *gradientCqDotDot_matrix_3; // C
	double alpha; // baumgarte stabilization coefficient
	double beta; // baumgarte stabilization coefficient
	double *right_matrix;
	gsl_vector *x;

	// ring variables
	float inner_radius; // minor radius
	float outter_radius; // major radius
	int num_major;
	int num_minor;
	point ring_position;

	//void ParticleSystemInit();
	// compute constraint
	double LinkConstraint(int index);
	double RingConstraint(int index); // just for the last particle
	// compute gradient of the constraint
	void LinkGradient(int index); 
	void RingGradient(int index);
	// compute dC' / dq
	void ComputeLinkCDot(int index);
	void ComputeRingCDot(int index);
	// mass matrix M
	void MassMatrix();
	// gradient C matrix
	void GradientCMatrix();
	// gradient C' matrix
	void GradientCDotMatrix();
	// transposed gradient C matrix
	void TransposedGradientCMatrix();
	// combine all four matrix to left matrix
	void LeftMatrix();
	// external force matrix
	void ExternalForceMatrix();
	// matrix multiplication
	void MatrixMultiplication(double *matrix_top, double *matrix_new, double *matrix_out);
	// build right matrix
	void RightMatrix();

	// render function
	// draw particle as a sphere
	void glDrawSphere(float s_radius, int s_subdivisions);
	// draw links as cylinders
	void glDrawCylinder(float c_radius, int c_subdivisions);
	// draw the ring
	void glDrawRing();
};
#endif