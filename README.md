# SPH CODE

Build the project using the standard CMake workflow (please build in release mode):

```
mkdir build
cd build
cmake ..
make
./fluid_project
```

# Dependencies
This project needs the libigl and polyscope libraries.
 * libigl: https://github.com/libigl/libigl
 * polyscope: https://polyscope.run/
 * OpenMP

Libigl should get installed automatically. You will need to download the Polyscope library yourself, and set the `POLYSCOPE_DIR` environment variable to the path where you installed Polyscope.

# Contents

The main contents of the simulation and demo are located in the file `main.cpp`. This file contains the code that will simulate the dam breaking demo. This code is a heavily modified version of the original fluid simulation code. It still uses the same code syntax and layout. We have a tight loop which we call the necessary functions documented in our implementation. The state of simulation is held in a simdata struct which contains all the necessary eigen vectors and matrices of the simulation. There are several function used to perform initialization of the particles, and polyscope which are modified from the original fluid project. 

The main function implements the polyscope callback. It first sets up the necessary constants and calls makegrid which generates the particles at their appropiate positions for the start of the simulation. We then use a polyscope pointcloud to represent the particles. Additionally, we have code that will detect the position of the mouse and set a global variable for mouseClicked to true to allow use to implement a mouse force in the simulation. The main loop also allows the user to tweak the constants stored in simdata. Some of the constants and changes will reset the simulation. Then within the polyscope call back we execute one step of the simulation.

The integrate function will apply velocity verlet time integration to the particles. Additionally, it will add the external forces of gravity, the mouse attraction force, and a velcotiy based impulse to keeep particles within the box of our simulation.

The compute density function will loop over all particles and their neighbors to calculate density based on our weighting function. This is optimized to perform in parallel using omp. Additionally, we try to forgoe computing square roots until we identify the particle as a neighbor futher reduces time complexity.

The compute pressure function will calculate the necessary pressure values in parallel using the gas equation.

The pressure force function calculate the force of pressure by iterating over all particles and their neighbors and using the equation for the force of pressure to sum up a pressure force on each particle.

The viscosity force similar to pressure force will calculate a force for viscosity using the particles and their neighbors.

# Extra Credit Honor Statement

On our honor, Nalin Mahajan and Vineeth Bandi have submitted the online course instructor survey.
