/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeAcceleration(struct world * jello, struct point a[8][8][8]);
void hookForce(const point& a, const point& b, double kHook, point& c);
void hookForceDiagonal2D(const point& a, const point& b, double kHook, point& f);
void hookForceDiagonal3D(const point& a, const point& b, double kHook, point& f);
void hookForceBend(const point& a, const point& b, double kHook, point& f);
void hookForceCollision(const point& a, const point& b, double kHook, point& n, point& f);
void dampingForce(const point& a, const point& b, const point& va, const point& vb, double kDamp, point& c);
void dampingForceCollision(const point& a, const point& b, const point& va, const point& vb, double kDamp, point& n, point& f);
void structualForce(world *jello, int i, int j, int k, point& a);
void shearForce(world *jello, int i, int j, int k, point& a);
void bendForce(world *jello, int i, int j, int k, point& a);
void internalForce(world *jello, int i, int j, int k, point& a);
void externalForce(world *jello, int i, int j, int k, point& a);
void collisionForce(world *jello, int i, int j, int k, point& a);
void inclinedPlaneSide(world *jello, bool& sameSide, point& normal);
void planeCollisionForce(world * jello, int i, int j, int k, const bool& s, const point& pn, point &a);
// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);

#endif

