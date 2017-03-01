# argonLJMD

Usage: bash ./pipeline.sh
	 	 	 	
## Introduction

Molecular dynamics (MD) involves simulating the physical movements and interactions of atoms and molecules in a computer. Such simulations allow researchers to understand properties of a system of molecules or atoms, from first principles. A fundamental requirement of MD simulations is provided by the definition of the interaction potential. One can describe pairwise interatomic forces with Lennard-Jones potential, for noble gases such as Argon whose atoms behave approximately like hard spheres that attract one another with Wan der Waals forces:
	
LJ (r) =V(r)= 4[(r)12-(r)6]

Where r is the distance between two atoms’ centers ε is the energy-parameter, the strength of the potential energy, and σ is the atom’s size-parameter. 

The primary goal of this study is to approximate macroscopic properties of a system of gaseous Argon atoms: (i) diffusion coefficient, and (ii) temperature. In this study, measurements of these properties were achieved by implementing an MD simulation, for argon atoms interacting according to the Lennard-Jones potential, and using the velocity verlet integrator to solve the equations of motion. The molecular dynamics of a system of 125, N, gaseous Argon atoms were simulated by implementing and program compiled in the C++ computing language (md.cc), plotting of data was implemented in Python (plot.py) and a BASH shell script to coordinate execution of md.cc, plot.py and to prepare the environent for multi-thread parallelization (pipeline.sh). 

## Reduced units

The theorem of corresponding states allows us to define a set of dimensionless reduced units. One particular implication of this theorem is that we can perform all simulations with ε=1.0 and σ=1.0, the so called Lennard-Jones units. Reducing values to dimensionless Lennard-Jones units simplifies calculations, and keeps values within a sensible dynamic range during the MD simulation, for example, the Lennard-Jones force between two particles, i and i is given by:
F i,j LJ= -Fj,iLJ=48r14-24r8

where r = |ri − rj|, and i and j correspond to pairs of particles. Henceforth, LJsignifies that a value is in Lennard-Jones units to be converted to absolute values for comparison with experimental values. 
Velocity verlet solver

The velocity verlet integrator was used to solve the equations of motion, to obtain values for positions, r, and velocities, through time, t:

r(t+t)=r(t)+v(t)t+12a(t)t+12a(t)t2

v(t+t)=v(t)+12[a(t)+a(t+t)]t

Velocity verlet integrates with low computational cost. Errors in this algorithm are of O(dt4), stable in MD applications, and in particular it rather successfully conserves energy of the system. 

## Periodic boundary conditions and volume

To avoid unusual interactions with a physical boundary, periodic boundary conditions were implemented, in which a particle on, say, the right boundary of the system interacts with the particles close to the left boundary of the system as though there would be a copy of the whole system (a periodic image of the system). This fixes the volume, V, of the system at L3, where L=50 is the length of a lattice edge V=503. 

Density (N/V) was kept at a constant value of 10-3. As such, the initialization phase of the MD simulation, 125 argon atoms are evenly spaced within the lattice, as shown in Figure.2. 

Figure.1. Evenly spaced Argon atoms initialized for MD simulation.

Diffusion coefficient

One way to obtain a self-diffusion coefficient, D, is from position data of particles, which our MD simulation provides. D is proportional to the mean square displacement, MSD of a particle per time, t, where t approaches infinity, and where the MSD is the ensemble average over all particles and all time origins, t0 , is given by:

DLJ = 16 t<r(t0+t)-r(t0)2>t

Fig.2. Reveals how DLJevolves up to large values of t, simulated for dt=0.001 and 5 million, k, time-steps:

Figure.2. DLJthrough time for four independent MD simulations with different random initial velocities, with final values (in LJ units) of 316.93 (red), 291.58 (green) and 260.78 (navy blue), 275.19 (cyan).
Temperature

If the system is in thermal equilibrium, then Boltzmann’s Equipartition Theorem relates the temperature to kinetic energy:
TLJ=13(N-1)i=1Nvi2 
The MD simulation initiates velocities randomly, with a maximum of Vmax=2.7 and all subsequent change in velocities are the result of the Lennard-Jones forces between particles. Sometimes collisions may occur, which lead to sharp drops in pairwise velocity and thus a lower temperature of the system, but overall temperature remained fairly consistent at around 2.43 LJ units.

Figure.3. Minor temperature fluctuations around ~2.43 (LJ units)

## Algorithm performance and optimisation

Average runtime for 5 million timesteps and 125 molecules is 25 minutes (intel i3 processor).

Serial speed-up was achieved by targeting pairwise particle force calculations. Computing these forces according to Lennard-Jones potentials is the most time-expensive component of the MD program. Compute time for calculating pairwise forces between n particles is O(n2). These calculations were largely circumvented by skipping all interactions greater than a distance threshold (rcut), thereby vastly reducing the number of force calculations. We chose an rcut threshold of 3 (three times the size parameter of the atom).

Memory conservation was achieved by ensuring time-dependent value calculations were only stored into memory for the current timestep and the previous time-step, writing all values to output files, changing memory usage from O(nk) to O(k). A downside of this technique is that greater runtimes occur due to I/O file operations that take place at every t.

Parallelization with OpenMP was implemented on the pairwise Lennard-Jones force calculations. Splitting the n2 pairs of calculations amongst three threads, on a two core (four thread) machine. Three threads yielded optimal runtimes saving on average 26% of the runtime. Using either a single thread or all four threads led to similar runtimes, both worse than using three.  

Goal: Use Molecular Dynamics to simulate simple system of argon atoms 
Supervision: Numerical Methods towards Exascale Workshop Wuppertal
