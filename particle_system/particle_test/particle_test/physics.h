/*

  USC/Viterbi/Computer Science
  "Particle System" Assignment 3 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_
#include "particle.h"

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(ParticleSystem *ps);
void RK4(ParticleSystem *ps);

#endif

