# TURBINS: A code to simulate fluid flows in oceans, atmosphere, and many other applications

This is the source code (written in C) of my code that I developed during my PhD studies at the Mechanical Engineering Department, University of California, Santa Barbara. 

TURBINS (**Turb**idity currents vis **I**mmersed-boundary **N**avier-**S**tokes) is capable of simulating gravity currents, particle-laden flows, sediment transport, and similar flow flows in geophysical fluid flows. The details of the code, validation, and parallelism are provided in this [scientific paper](http://dx.doi.org/10.1016/j.compfluid.2010.11.023). 

TURBINS is a 3D, parallel code which uses sophisticated domain decompostion method, MPI and distributed arrays for parallelism. Iterative linear solvers are preconditioned to solve the discretized momentum and concentration transport equations. User can set different solvers and preconditioners in Solver.c considering the problem at hand. 

TURBINS has more than ~25,000 lines of C-code with an object-oriented paradigm. The post-processing code (not available in this repository) is written in C++ and converts the raw binary output files to [VTK format](http://www.vtk.org/).

I have utilized TURBINS to run on O(1000) processors on different supercompting facilities and gained good speed-up. I had to transfer, polish, extract statistical quantities from O(200) Tera-bytes of resulting data. 

TURBINS has been utilized in more than 15 journal (peer-reviewed) publications over the past few years by scientists worldwide. For a complete list of punlications, please visit the [CFD-Lab](https://sites.google.com/site/ucsbcfdlab/) at UCSB. 


**Requirements:**

- [petsc](https://www.mcs.anl.gov/petsc/) (version 3.1) 
- [hypre](http://acts.nersc.gov/hypre/) (version 2.6 and above)
- suitable [MPI package](https://www.open-mpi.org/)

**Visualization:**
- Any software supporting VTK format, e.g. [ParaView](https://www.paraview.org) (or [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit))
- Rendering of some of the shown movies in the [gallery](/gallery) sectio is performed by ParaView and Python. 

Note: Since this code is a research code, the main.c file, Makefile, and input files are not uploaded here. If you are interested in using this code, please contact me personally via mmnasr@gmail.com

My personal interests are: C/C++, python, scientific computing, parallel programming, data science, machine learning
