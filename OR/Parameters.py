"""
Plug in parameters
"""
import Constant_JGRpaper as c
import matplotlib.pyplot as plt

#Solver parameters
t0= 0 #starting time
tf= 50 *c.secinday #(s) indicate number of days until solver stops
h0= 795 #(m) initial head
S0= 1.2 #(m) initial channel cross-section area

R_mean= 3 # m^3/s, mean discharge into moulin
R_min= 2.6 #Minimum daily recharge (used for sine func)

R_period= c.secinday #Period for oscillatory Recharge
H= 1000 # m, thickness of the ice
L= 30000 # m, length of the conduit
Pi= c.rhoi*c.g*H #Units?? Ice pressure 
#r_base = 5
#p['r_hmiddle']= 5




'''Figure parameters'''

fontsize_label = 8
fontsize_tick = 8

dpi=300

