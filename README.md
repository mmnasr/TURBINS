# TURBINS: A code to simulation fluid flows in ocean and atmosphere

This is the source code (written in C) of my code that I developed during my PhD studies at the Mechanical Engineering Department, University of California, Santa Barbara. 

TURBINS (**Turb**idity currents vis **I**mmersed-boundary **N**avier-**S**tokes) is capable of solving gravity currents, particle-laden flows, and similar flow fields in geophysical fluid flows. The details of the code, validation, and parallelism are provided in this [scientific paper](http://dx.doi.org/10.1016/j.compfluid.2010.11.023). 

TURBINS has more than ~25,000 lines of C-code with an object-oriented paradigm. Post-processor code (not available in this repository) is C++ code which converts raw binary output files to  [VTK format](http://www.vtk.org/).

TURBINS is a 3D, parallel code which uses sophisticated domain decompostion method, MPI and distributed arrays for parallelism. Iterative linear solvers are preconditioned to solve the discretized momentum and concentration transport equations. User can set different solvers and preconditioners in Solver.c considering the problem at hand. 

TURBINS has been utilized in more than 15 journal (peer-reviewed) publications over the past few years by scientists worldwide. For a complete list of punlications, please visit the [CFD-Lab](https://sites.google.com/site/ucsbcfdlab/) at UCSB. 


**Requirements:**

- [petsc](https://www.mcs.anl.gov/petsc/) (version 3.1) 
- [hypre](http://acts.nersc.gov/hypre/) (version 2.6 and above)
- suitable [MPI package](https://www.open-mpi.org/)

**Visualization:**
- Any software supporting VTK format, e.g. [ParaView](www.paraview.org) (or [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit))
- Rendering of some of the shown movies in the [gallery](/gallery) sectio is performed by ParaView and Python. 

Note: Since this code is a research code, the main.c file, Makefile, and input files are not uploaded here. If you are interested in using this code, please contact me personally via mmnasr@gmail.com

Interests: C/C++, python, scientific computing, parallel programming, data science, machine learning
