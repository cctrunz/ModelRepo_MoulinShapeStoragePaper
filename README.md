# Model Repository
 This repository is a supplement to the paper "Modeling the impact of moulin shape on englacial hydrology", currently submitted to JGR earth processes

The folders called Constant_Qin and Oscillating_Qin hold the code used to simulate the head in moulins with a constant meltwater input and a oscillating meltwater input respectively.

This repository contains all the functions and ipython notebook to create the figures for the paper

The model runs with a standard conda environment. In case it does not run, you can use this environment:
> conda create --name <env> --file package-list.txt


## Simulations with constant input, Section 3.1
Observe effect of moulin shape on moulin head oscillation towards equilibrium when the forcing changes abruptly.

### [Code to reproduce the subpanles in Figure 2](Constant_Qin/Figure2.ipynb) -- Simulations for 4 different shapes, with the same forcing and same boundery conditions

### [Code to reproduce Figure 3](Constant_Qin/Figure3.ipynb) -- Simulations for the same 4 different shapes as Figure 2, but for several position along an idealized ice sheet

## Simulations with sinusoidal input, Section 3.2 and 3.3 

### [Code to reproduce the subpanels in Figure 4](Constant_Qin/Figure3.ipynb)
