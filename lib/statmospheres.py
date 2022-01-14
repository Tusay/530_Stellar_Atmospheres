# This code will be the function repository for useful functions from 530 Stellar Atmospheres.
# Include it in python with the import module feature (use sys to pull it from another directory if need be)

from astropy import units as u
from astropy import constants as const
import numpy as np

# Calculate the Blackbody Planck Function for one frequency (nu) and one temp (T). 
# Optional coefficient A. 
# cgs units.
@u.quantity_input(T=u.K,nu=u.Hz)
def one_Planck(T,nu, A=1) -> (u.erg * u.cm**-2 * u.Hz**-1 * u.s**-1 * u.sr**-1):
    B_nu = A * (2 * const.h.cgs * nu**3 / (const.c.cgs)**2) / ((np.exp(const.h.cgs*nu/(const.k_B.cgs*T)))-1) / u.sr
    return B_nu

# Calculate the Blackbody Planck Function for an array of frequencies (nu) 
# and one temp (T) using the function one_Planck() defined above. 
# cgs units
def Planck_array_at_temp(T,nu_array,A=1) -> (u.erg * u.cm**-2 * u.Hz**-1 * u.s**-1 * u.sr**-1):
    B_nus=[]
    for nu in nu_array:
        B_nus.append(one_Planck(T,nu,A).value)
    return B_nus*(u.erg * u.cm**-2 * u.Hz**-1 * u.s**-1 * u.sr**-1)