"""
To get a age grid
"""

import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from scipy.stats import skewnorm, exponnorm

def CranmerSaar2011_eqn36(Teff):
    """
    https://ui.adsabs.harvard.edu/abs/2011ApJ...741...54C
    """
    tau = 314.24 * np.exp(-1*(Teff / 1952.5) - ((Teff/6250)**18)) + 0.002
    return tau
def bprp2teff(x):
    Teff = -4.16584*10e1\
    +3.97800*10e3*x\
    -8.41905*10e3*x**2\
    +8.52039*10e3*x**3\
    -4.82259*10e3*x**4\
    +1.55985*10e3*x**5\
    -2.69476*10e2*x**6\
    +1.92865*10e1*x**7
    return Teff
def age_grids(x, c0,c1,c2,c3,c4,c5,c6, Ro1, Ro2, skw, samp):
    bprp, age = x
    log_bprp = np.log10(bprp)
    log_age = np.log10(age*10**6)
    Teff =  bprp2teff(bprp)
    tau = CranmerSaar2011_eqn36(Teff)
    P1 = Ro1 * tau
    P2 = Ro2 * tau
    t1 = 10**(np.log10(P1)/c6-(c0+c1*log_bprp+c2*log_bprp**2+c3*log_bprp**3+c4*log_bprp**4+c5*log_bprp**5)/c6)
    t2 = 10**(np.log10(P2)/c6-(c0+c1*log_bprp+c2*log_bprp**2+c3*log_bprp**3+c4*log_bprp**4+c5*log_bprp**5)/c6)
    dP = skewnorm(skw, loc=(t1+t2)/2, scale=(t2 - t1)/2).pdf(age*10**6)
#     plt.hist(dP)
    dP = dP*(t2 - t1)/2 * np.sqrt(2*np.pi) / 2 * np.exp(1/(2*np.abs(skw))) * ((P1-P2)*samp)
# 
    P = dP + 10**((log_age *c6) + (c0+c1*log_bprp+c2*log_bprp**2+c3*log_bprp**3+c4*log_bprp**4+c5*log_bprp**5))
    return P
def get_age_grids(age_range=np.linspace(700, 4000, 50)):
    p = [-3.09288236e+00,  3.94419542e-01, -3.13872184e+00,  1.62308904e+01,
        -8.19160418e+00, -1.73019171e+01,  4.69446261e-01,  6.09199538e-01,
         9.10678414e-01, -7.04677735e-02,  1.08314195e-03]
    bprp = np.linspace(0.5, 2.5, 1000)
    ini_grid = pd.DataFrame(data=None,columns=['bp_rp0','age','Prot'])
    for i in range(len(age_range)):
        dic = {'bp_rp0': bprp, 'age': age_range[i],'Prot':age_grids((bprp, np.ones(len(bprp))*age_range[i]), *p)}
        grid = pd.DataFrame(dic)
        ini_grid = pd.concat([ini_grid,grid],axis=0,ignore_index=True)
    return ini_grid

