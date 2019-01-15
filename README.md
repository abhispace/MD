# MD
Molecular Dynamics Simulation

Introduction
The objective of this project is to study the extensivity of internal energy for a dusty plasma system of N particles that interact through a Yukawa potential. We analyze the evolution of potential and total energy per particle with time using the Verlet Algorithm. We apply different values of screening parameter that modifies the strength of potentials experienced by the particles. We show that internal energy is an extensive property for high values of the screening parameter or conversely, we show that Yukawa potential is a short range potential for high values of the screening parameter.

Codes:
rand_pos.f90: creates the initial positions and velocities
main.f90:Reads the input files
sub_gen_particles.f90: reads input positions and velocities
force.f90: calculates the force on each particle
integrator.f90: uses the Verlet algorithm to calculate new positions

Instructions:
Run rand_pos.f90 first to get the initial positions and velocities of particles. Change the number of particles according to your needs.

Then in terminal run the rest of the four forrtan files together.
