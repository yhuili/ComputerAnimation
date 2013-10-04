ComputerAnimation
=================

This is a USC CSCI520 HW Repo
1. Jello Cube
In this project, we need to simulate a jello cube (8*8*8) that moves around within a bounding box under an arbitrary non-uniform force field.
In order to simulate this jello cube, we use a 3D mass-spring network, which includes structural, shear and bend springs.
In addition, we also create an external force field to affect the status of the cube.
Finally, we need to simulate the collision detection and response for the jello cube hitting any of the six walls of the bounding box.
With the help of create world file, we can define the cube environment ourselves. Then, we read the world using input.cpp.


2. Mocap Player Motion Interpolation
In this project, we implement three interpolation schemes to interpolate human motion data obtained from an optical mocap system, 
including Bezier interpolation for Euler angles, SLERP and Bezier interpolation for quaternion.

3. Particle System
In this project, we use particle system to model a simple chain, with an emphasis on "hard" constraints.

Please check detailed info from the README file in each folder.
