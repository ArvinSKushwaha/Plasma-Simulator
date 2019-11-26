# Plasma-Simulator
This is code for a plasma simulator I am working on for tokamaks. It has a long way to go!

* math_obj.hpp: Contains useful mathematical objects
    * Vec3D: A vector class
    * SizeTuple: An unsigned integer class (to hold matrix sizes)
    * ScalarField: A 3D matrix (field) of scalars; it is a discretized representation of a scalar field function
    * VectorField: Tmamaxwhe Vector analog of the ScalarField. Represents a discretized representation of a vector field function.

* pure_field/main.cpp: Code for a fluid model simulation of plasma. Assumes non-relativistic situations.
    * Contains field initialization
        * NEEDS WORK (magneticField doesn't initialize properly)
    * Contains field propagation
        * Unsure if propagation works, more testing necessary

* particle_cell/main.cpp: Mostly empty code for PIC model simulation of plasma
    * Plan on implementing PIC model

* particle_cell/plasma_sim.hpp: Contains supporting functions and classes for particle_cell/main.cpp.
    * Needs a lot more work to finish functions
    * gammaFromV: Calculates lorentz factor for a velocity vector.
    * YeeLattice: Contains information on the magnetic and electric field in the simulation
    * Particle: A particle object for superparticle simulation.