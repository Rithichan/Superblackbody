# Superblackbody
Simulation for thermal radiation

## To do
Physics
1. Specifics on thermal fluctuations leading to radiation
2. Calibration of the C correction constant
3. Specifics on the derivation of K and other equations that may be relevant
Figure out how to simulate thermal radiation by adding thermal fluctuations in tidy3D
1. Currently have a method to generate random current density J(t) by taking random thermal fluctuations K(t), need to test this.
2. Find a method to create this custom source
Explicit method directly using Maxwell equations solving for Differential equation across the grid
1. Need to implement methods for calculating relevant quantities, for now have Poynting vector and Power flux (check definition).
2. Implement boundary conditions
3. More testing
4. Optimisation
