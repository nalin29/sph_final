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


