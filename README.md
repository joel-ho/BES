# BES

Basic Euler Solver, a lightweight, multi-threaded, three-dimensional unstructured finite volume code to solve the Euler equations. More information at the [project page](https://joelcfd.com/projects/bes).

BES uses [OpenMP](http://www.openmp.org/), [Eigen](http://eigen.tuxfamily.org/), [AMGCL](http://amgcl.readthedocs.io/en/latest/) and [Boost](https://www.boost.org/). It is released under version 3 of the GNU General Public License.

## Getting Started

Before deciding to compile the code, you can check out the binaries for [Windows 64 bit](https://drive.google.com/file/d/18ypL0WMQS3lU-KcvVytyWhVVcTkJy6Wc/view?usp=sharing) or [Linux 64 bit]().

If you want to compile, you just need to have the source file for Boost somewhere on your machine. Some header-only libraries from Boost is used in AMGCL, hence the dependency. The required Eigen and AMGCL libraries are provided in BES include directory, so no action is needed for that.

### Prerequisites

You can download Boost [here](https://www.boost.org/users/history/). Boost v1.66.0 was used in the development of BES.

MinGW-w64 v4.3.3 with gcc v7.2.0 was used in the development for windows.

You might also need [Gmsh](http://gmsh.info/) for meshing the test cases and your own cases in future. The results can be viewed with [Paraview](https://www.paraview.org/).

### Installing

The compilation of BES is from a simple `makefile`, Automake or Autoconf was not used. As such, you will need to manually input some options into the `makefile`.

First, enter your OS in line 2 of the `makefile`. For Windows:

```
OS = WINDOWS
```

For Linux:

```
OS = LINUX
```

Second, enter the path to Boost in your machine in line 36 of the `makefile`:

```
BOOST = path/to/boost
```

We are good to go! For Windows, we can compile in parallel with N threads by:

```
mingw32-make -j N
```

For Linux:

```
make -j N
```

The binary will be created in the `bin` folder. You can add it to your path if you desire.

## Running BES

We will go through an example with a two-dimensional NACA0012 test case. Starting from the BES folder, we navigate to the case folder:

```
cd test\twoDimNacaAirfoil
```

The Gmsh `.geo` files are provided in the test cases. You will need to create the mesh youself, i.e.

```
gmsh nacaAirfoil.geo -2
```

Then you can configure the case in the `.besi` file. For example, if you want to change the freestream Mach number, you can change the FREESTREAM_MACH option in `nacaAirfoil.besi`:

```
FREESTREAM_MACH = 0.5,0.05
```

This will set the freestream Mach to be 0.5 in the x-direction and 0.05 in the y-direction. Also, ensure that the mesh file is correctly named in the MESH_FILE option. The results will be written in the folder named in the OUTPUT_FOLDER option. The folder must be created before running the case or else BES will not output the results. Remember to save, and we can run BES now. We give our input file as the one and only argument to the program:

```
BES nacaAirfoil.besi
```

Now sit back and relax! Once done, you can view the results in the output folder with Paraview or any other viewer that can open `.vtk` files.

## Authors

* **Ho Mun Onn, Joel** - [JoelCFD](https://joelcfd.com/about)

## License

This project is licensed under version 3 of the GNU General Public License - see the [LICENSE](LICENSE) file for details.
