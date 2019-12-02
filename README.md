# Plasma-Simulator
This is code for a plasma simulator I am working on for tokamaks. It has a long way to go!
This repository also has code for an AI model that contains an RCNN (Recurrent Convolutional Neural Network) for simulations (Tensorflow backend).

* math_obj.hpp: Contains useful mathematical objects
    * Vec3D: A vector class in 3D space
    * Vec2D: A vector class in 2D space
    * SizeTuple: An unsigned integer class (to hold matrix sizes)
    * ScalarField: A 3D matrix (field) of scalars; it is a discretized representation of a scalar field function
    * VectorField: The Vector analog of the ScalarField. Represents a discretized representation of a vector field function.
    * ScalarPlane: A 2D matrix (field) of scalars; Is the 2D analog of the ScalarField
    * VectorPlane: A 2D matrix (field) of vectors; Is the Vector analog of the ScalarPlane

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

* plasma_AI/main.py: Contains a C-RNN for AI-based plasma simulations.
    * Needs to be trained on plasma simulation, but have to wait until the rest of the code is prepared.