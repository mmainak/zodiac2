# zodiac2
Description of the Code:  Zodiac is a pseudo spectral code with adjustable choices of domain periodicity. It is mainly developed for analyzing boundary layer turbulence with complex geometry. It uses a Rk3 solver, ADI for the the nonlinear terms and Multi-grid for convergence. It is MPI parallelised and the domain-decomposition  can be customised based on the needs. The structure of the code is briefly described here.  zodiac.F90: It is the  main program that initialise the code and get input parameters (input.dat, grid_def), calls in the solver, manages  output.  duct.F90: solves the RK3 steps of the NS equations, uses FFT, calls in ADI for the implicit non-linear part, uses Multi-grid for faster convergence, deals with different forcing (written inside the same files as subroutines, viz. wave_forcing)  boundary.F90: Provides the boundary conditions, in the present case it  has been used to study the melting problem in the ice-water interface. It can be customised based on the problem).  flow_statistics.F90: Deals with the tke and energy budget of the flow. It also stitches the decomposed domain from other processors.  flow_output.F90: Deals with the output, presently it can output in .plt (for tecplot, that need tecio64.a files), .vtk(for preview) and also for netcdf(to be added as a patch).  GRID_XY.dat:  It is the main grid file for the present problem, that can be customized.

Required Libraries:
gcc/4.9.3 (or higher)
netcdf-fortran
fftw2 (the code is not compatible with fftw3)
openmpi(1.10.2 or higher, for multiprocessing)
tecplot library file (tecio64.a)

The uses needs to create the blank directories: plane_data, plane_tke, last_saved, 
