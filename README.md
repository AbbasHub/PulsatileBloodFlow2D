# PulsatileBloodFlow2D
This is a Lattice-Boltzmann solver that was first written in Fortran 90 and then (very crudely) re-written in C++. It is an LBM code for simulation of pulsatile flow through compliant vessels.

The array shapes are kept as they were in the Fortran code without converting them to uni-vectors of pointers, which are ideal for an optimized C++ code. The C++ version of the code strives for readability rather than efficiency.
