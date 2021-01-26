import numpy as np

"""
Fixed Parameters
"""
rhow = 1000 #water density, kg/m3
rhoi =  910 #ice density $, kg/m3
g = 9.8 #Gravity
f = 0.1 #Darcy weissbach friction factor #different in Schoof(2010)
Lf = 3.32e5 #Latent heat of fusion (from schoof2010) J/kg  #(CT) 3.35e5 in schoof(2010)
A = 6e-24 #Pa^-3 s^-1
n = 3 #flow law parameter

#Calculated parameters
C1 = 1/(rhoi*Lf)
C2 = 2*A*n**-n
C3 = ( 2.**(5./4.)/np.pi**(1./4.) * np.sqrt(np.pi/(np.pi + 2.)) )/ np.sqrt(rhow*f) #(2**(1/4) * (np.pi+2)**(1/2)) / (np.pi**(1/4) * (rhow*f)**(1/2))
secinday = 24*3600