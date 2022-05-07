# Fluid simulation starter code

Build the project using the standard CMake workflow (please build in release mode):

```
mkdir build
cd build
cmake ..
make
```

# Dependencies
This project needs the libigl and polyscope libraries.
 * libigl: https://github.com/libigl/libigl
 * polyscope: https://polyscope.run/
 * OpenMP

Libigl should get installed automatically. You will need to download the Polyscope library yourself, and set the `POLYSCOPE_DIR` environment variable to the path where you installed Polyscope.

# Contents

The main contents of the simulation and demo are located in the file `main.cpp`. This file contains the code that will simulate the dam breaking demo. This code is a heavily modified version of the original fluid simulation code. It still uses the same code syntax and layout. We have a tight loop which we call the necessary functions documented in our implementation. The state of simulation is held in a simdata struct which contains all the necessary eigen vectors and matrices of the simulation. There are several function used to perform initialization of the particles, and polyscope which are modified from the original fluid project. 


