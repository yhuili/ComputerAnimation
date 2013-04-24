<Please submit this file with your solution.>

CSCI 520, Assignment 1

Name: Hui Li
ID: 3194819503
Email Address: hli997@usc.edu
VS Version: Microsoft Visual Studio 2010 Professional Version 10.0.40219.1 SP1Rel

================

<Description of what you have accomplished>
Objectives:
1. In this project, we need to simulate a jello cube (8*8*8) that moves around within a bounding box under an arbitrary non-uniform force field.
2. In order to simulate this jello cube, we use a 3D mass-spring network, which includes structural, shear and bend springs.
3. In addition, we also create an external force field to affect the status of the cube.
4. Finally, we need to simulate the collision detection and response for the jello cube hitting any of the six walls of the bounding box.
5. With the help of create world file, we can define the cube environment ourselves. Then, we read the world using input.cpp.


Some Instructions:
0. If you want to compile it under vs2010, first, please configue GLUT files, add glut32.lib opengl32.lib to Additional Dependencies under Configuration Properties/Linker/Input. Second, add GLUT files directory to Include Directories and Library Directories under Configuration Properties/VC++ Directories
1. Please add world file such as world\jello.w under Configuration Properties/Debugging/Command Arguments.
2. Please change the solution configuration to Release, and run it.

What I have done:
1. There are several additional functions to help finish the computation, they are:
// calculate the hook force between point a and its neighbour point b
void hookForce(const point& a, const point& b, double kHook, point& c);

// calculate the hook force between point a and its diagonal neighbour point b on surface
void hookForceDiagonal2D(const point& a, const point& b, double kHook, point& f);

// calculate the hook force between point a and its diagonal neighbour point b in cube
void hookForceDiagonal3D(const point& a, const point& b, double kHook, point& f);

// calculate the hook force between point a and its the next neighbour point b
void hookForceBend(const point& a, const point& b, double kHook, point& f);

// calculate the hook force for collision springs
void hookForceCollision(const point& a, const point& b, double kHook, point& n, point& f);

// calculate the damping force between point a and point b
void dampingForce(const point& a, const point& b, const point& va, const point& vb, double kDamp, point& c);

// calculate the damping force for collision springs
void dampingForceCollision(const point& a, const point& b, const point& va, const point& vb, double kDamp, point& n, point& f);

// calculate the structual force between point p at (i, j, k) and its neighbour point b
void structualForce(world *jello, int i, int j, int k, point& a);

// calculate the shear force between point p at (i, j, k) and its diagonal neighbour point b
void shearForce(world *jello, int i, int j, int k, point& a);

// calculate the bend force between point p at (i, j, k) and its following neighbour point b
void bendForce(world *jello, int i, int j, int k, point& a);

// calculate the internal force, including structual, shear, and bend force, add them to point p at (i, j, k)
void internalForce(world *jello, int i, int j, int k, point& a);

// calculate the external force field, add it to point p at (i, j, k) 
void externalForce(world *jello, int i, int j, int k, point& a);

// calculate the collision response force by the bounding box walls
void collisionForce(world *jello, int i, int j, int k, point& a);

2. I also changed the doIdle function, capturing the screen every 5 frames.

<Also, explain any extra credit that you have implemented.>
For extra credit, I drew an inclined plane in the bounding box, and add a collision detection and response force for it.
I first check if the jello cube and the plane normal are on the same side using function:
void inclinedPlaneSide(world *jello, point& p, bool& sameSide, point& normal)

0. if they are not on the same side, I changed the direction of the normal.
1. if they are on the same side, I go to the loop, check the first time that the cube point and the plane normal are on different sides.
2. then, I calculate the response force using function:
void planeCollisionForce(world * jello, int i, int j, int k, const bool& s, const point& pn, point &a)