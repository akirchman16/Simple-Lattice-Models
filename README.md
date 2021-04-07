# Simple Lattice Models
 Original lattice models used for learning basics of the model.

This folder contains MATLAB codes which utilize a simple 1-dimensional Monte Carlo model to describe the ssDNA
saturation as RAD51 binds to and unbinds from the ssDNA. The three primary codes in this folder include:

	1. NonInteracting_LatticeModel_Scatchard
	2. NonInteracting_LatticeModel
	3. Interacting_LatticeModel_Scatchard

All of the above codes utilize the same process:
	This code uses probabilities of binding and unbinding based on ODE models of chemical reactions. With each
	iteration of the model, each available binding location is checked, in a random order, and tested for 
	binding and then each bound protein is checked and tested for unbinding. An "available location" is
	classified as a location that is currently free and the next (n-1) locations on the DNA lattice are free,
	where n is the length of the protein. With each check at each location a random number, r, is selected
	between 0 and 1. At the same time, a probability or binding is determined using the following equations:
			P(binding) = k_on*L*dt
			P(unbinding) = k_off*dt
	where k_on and k_off are the kinetic association constants for binding and unbinding, respectively, L is
	the free protein concentration in solution, and dt is a small time step between each iteration. All of
	these values are input paramaters for the code. The random number, r, is compared to this probability and
	if it is smaller than the probability the corresponding reaction is carried out (a protein is bound at that
	location, or the chosen protein is unbound). After checking each location and protein on the DNA lattice,
	time is advanced by the small time step, dt. The process is then repeated at the new time value. Over time,
	the saturation level of the DNA lattice is recorded. This process repeats many times until an equilibrium
	is reached in the system.

The key differences between the 3 primary codes is whether they display pure results or whether they compare the
model the previous results, Scatchard plots, and whetehr they include cooperativity of RAD51 proteins.
