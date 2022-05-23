# Model Repository
This repository is a supplement to the paper "Modeling the influence of moulin shape on subglacial hydrology", currently submitted to JGR earth processes

Trunz, C., Covington, M. D., Poinar, K., Andrews, L. C., Mejia Zimmerman, J., & Gulley, J. (submitted 2022). Modeling the influence of moulin shape on subglacial hydrology. In review at JGR earth process.


The folders called Constant_Qin and Oscillating_Qin hold the code used to simulate the head in moulins with a constant meltwater input and a oscillating meltwater input respectively.

This repository contains all the functions and ipython notebook to create the figures for the paper

The model runs with a standard conda environment. In case it does not run, you can use this environment:
> conda create --name <env> --file package-list.txt


## Simulations with constant input, Section 3.1
Observe effect of moulin shape on moulin head oscillation towards equilibrium when the forcing changes abruptly.
The model and the plotting function are held in this python file: [FR_mainCode_plus.py](onstant_Qin/FR_mainCode_plus.py)

### [Code to reproduce the subpanels in Figure 2, S2 and S4](Constant_Qin/Plot_Figure2-S2-S4.ipynb) -- Simulations for 4 different shapes, with the same forcing and same boundery conditions

### [Code to reproduce Figure 3](Constant_Qin/Figure3.ipynb) -- Simulations for the same 4 different shapes as Figure 2, but for several position along an idealized ice sheet

## Simulations with sinusoidal input, Section 3.2 and 3.3 
For those simulations, the model coded as a function inside the notebook directly. 
Here you can find the [constants](Oscillating_Qin/Constant_JGRpaper.py) and the [parameters](Oscillating_Qin/Parameters.py) we used for the simulations.

### [Code to reproduce the subpanels in Figure 4](Oscillating_Qin/OscillationRecharge_Figures-Parameters-For-JGR.ipynb) -- Simulations for 5 different shapes.
  


## Simulations comparison between a 0D and a discretized 1D subglacial channel model, with an oscillating input. 
  The 1D model for the subglacial channel is [onedim_channel_model.py](Compare_0D_1D/onedim_channel_model.py)
### [Code to reproduce the supporting Figures S6 and S7](Compare_0D_1D/Supplemental_figures_S06_S07.ipynb) 
