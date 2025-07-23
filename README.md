# Superblackbody
Simulation for thermal radiation with tidy3D for wide and narrowband thermal emitter

# To do
## 1. Optimisation of credits
The simulation right now has a lot of redundant space below the metasurface by moving the structure lower into the simulation domain we can probably reduce credit costs. Right now one run of the simulation is around 0.6-0.8 credits. If possible, you may try to find other optimisations from your experience working with Fangrus code.
## 2. Combine both jupyter notebooks into one and remove redundancy
Both notebooks have a lot of things that are completely identical. Clean up the code and combine them into one.
## 3. Try other frequency widths, distributions and locations of the dipoles
Currently, the radiation is simulated by scattering random dipoles around the structure uniformly. The frequency of each dipole is distributed such that the frequencies match that of the blackbody curve with a large number of dipoles with random phases. There really isn't any analytical basis for this, which I am not too happy about, but it works. I am currently reading to find some backing for this methodology, but for now, you can play around with the frequency distributions and distributions of dipole locations.
Some of my suggestions:
- Try uniformly distributed frequencies for the dipoles (I have tried this and it works well too, but see if it works now)
- Try having all the dipoles have the same central frequency and a large frequency width.
- Try distributing the dipoles in another way (Joel and I have discussed this, maybe try having a certain number of dipoles in each layer or having the dipoles at the surface only)
Play around with different settings and see what works best. I have attached some references for the spectrum that you should get.
## 4. Run a large batch simulation
Once you have optimised the code to be cheaper in credits, you should run a large number of simulations, around a few hundred if possible. The process is inherently stochastic, so running more simulations and taking the sum of the flux should give us a more accurate result. Right now, with around 10 or so simulations, the main peak is clearly present, but the side peaks can vary abit due to randomness.
