/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include <math.h>

/* performs one step of Euler Integration */
/* as a result, updates the particle system structure */
// copy from HW1
void Euler(ParticleSystem *ps)
{
	ps->ComputeAcceleration();
	int num = ps->particle_num;
	for (int i = 1; i < num; i++)
	{
		// z is always 0
		ps->particle_list[i].p_position.x += ps->dt * ps->particle_list[i].p_velocity.x;
		ps->particle_list[i].p_position.y += ps->dt * ps->particle_list[i].p_velocity.y;
		ps->particle_list[i].p_velocity.x += ps->dt * ps->particle_list[i].p_acceleration.x;
		ps->particle_list[i].p_velocity.y += ps->dt * ps->particle_list[i].p_acceleration.y;
	}
}

/* performs one step of RK4 Integration */
/* as a result, updates the particle system structure */
// copy from HW1
void RK4(ParticleSystem *ps)
{
	int num = ps->particle_num;
	point *F1p = new point[num];
	point *F1v = new point[num];
	point *F2p = new point[num];
	point *F2v = new point[num];
	point *F3p = new point[num];
	point *F3v = new point[num];
	point *F4p = new point[num];
	point *F4v = new point[num];
	ParticleSystem *buffer = ps; // make a copy
	point *buffer_p = new point[num];
	point *buffer_v = new point[num];

	ps->ComputeAcceleration();
	for (int i = 1; i < num; i++)
	{
		//buffer->particle_list[i].p_acceleration;
		pMULTIPLY(ps->particle_list[i].p_velocity, ps->dt, F1p[i]);
		pMULTIPLY(ps->particle_list[i].p_acceleration, ps->dt, F1v[i]);
		pMULTIPLY(F1p[i], 0.5, buffer_p[i]);
		pMULTIPLY(F1v[i], 0.5, buffer_v[i]);
		pSUM(ps->particle_list[i].p_position, buffer_p[i], buffer_p[i]);
		pSUM(ps->particle_list[i].p_velocity, buffer_v[i], buffer_v[i]);
	}

	ps->ComputeAcceleration();
	for (int i = 1; i < num; i++)
	{
		// F2p = dt * buffer.v;
		pMULTIPLY(buffer_v[i], ps->dt, F2p[i]);
		// F2v = dt * a(buffer.p,buffer.v);     
		pMULTIPLY(ps->particle_list[i].p_acceleration, ps->dt, F2v[i]);
		pMULTIPLY(F2p[i], 0.5, buffer_p[i]);
		pMULTIPLY(F2v[i], 0.5, buffer_v[i]);
		pSUM(ps->particle_list[i].p_position, buffer_p[i], buffer_p[i]);
		pSUM(ps->particle_list[i].p_velocity, buffer_v[i], buffer_v[i]);
	}

	ps->ComputeAcceleration();
	for (int i = 1; i < num; i++)
	{
		// F3p = dt * buffer.v;
		pMULTIPLY(buffer_v[i], ps->dt, F3p[i]);
		// F3v = dt * a(buffer.p,buffer.v);     
		pMULTIPLY(ps->particle_list[i].p_acceleration, ps->dt, F3v[i]);
		pMULTIPLY(F3p[i], 0.5, buffer_p[i]);
		pMULTIPLY(F3v[i], 0.5, buffer_v[i]);
		pSUM(ps->particle_list[i].p_position, buffer_p[i], buffer_p[i]);
		pSUM(ps->particle_list[i].p_velocity, buffer_v[i], buffer_v[i]);
	}

	ps->ComputeAcceleration();
	for (int i = 1; i < num; i++)
	{
		// F4p = dt * buffer.v;
		pMULTIPLY(buffer_v[i], ps->dt, F4p[i]);
		// F4v = dt * a(buffer.p,buffer.v);     
		pMULTIPLY(ps->particle_list[i].p_acceleration, ps->dt, F4v[i]);

		pMULTIPLY(F2p[i], 2, buffer_p[i]);
		pMULTIPLY(F3p[i], 2, buffer_v[i]);
		pSUM(buffer_p[i], buffer_v[i], buffer_p[i]);
		pSUM(buffer_p[i], F1p[i], buffer_p[i]);
		pSUM(buffer_p[i], F4p[i], buffer_p[i]);
		pMULTIPLY(buffer_p[i], 1.0 / 6, buffer_p[i]);
		pSUM(buffer_p[i], ps->particle_list[i].p_position, ps->particle_list[i].p_position);

		pMULTIPLY(F2v[i], 2, buffer_p[i]);
		pMULTIPLY(F3v[i], 2, buffer_v[i]);
		pSUM(buffer_p[i], buffer_v[i], buffer_p[i]);
		pSUM(buffer_p[i], F1v[i], buffer_p[i]);
		pSUM(buffer_p[i], F4v[i], buffer_p[i]);
		pMULTIPLY(buffer_p[i], 1.0 / 6, buffer_p[i]);
		pSUM(buffer_p[i], ps->particle_list[i].p_velocity, ps->particle_list[i].p_velocity);
	}

	delete[] F1p;
	delete[] F1v;
	delete[] F2p;
	delete[] F2v;
	delete[] F3p;
	delete[] F3v;
	delete[] F4p;
	delete[] F4v;
	delete[] buffer_p;
	delete[] buffer_v;
	return;  
}
