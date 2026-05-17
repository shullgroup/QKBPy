# updated version of this file is maintained at
# https://github.com/shullgroup/QKBPy/blob/main/test/models.py

import numpy as np

# gaussian with baseline
def Gaussian(x, ctr, amp, wid, baseline=0):
        y = np.zeros_like(x)
        y = y + baseline + amp * np.exp(-((x - ctr)/(2 * wid))**2)
        return y

# arrhenius
def Arrhenius(T, A, Ea):

    return A * np.exp(Ea / (8.314 * T))

# vft
def VFT(T, A, B, Tinf):
      
      return A * np.exp(B / (T - Tinf)) 

# ln versions
def ln_Arrhenius(T, A, Ea):
      
      return np.log(A) + Ea / (8.314 * T)

def ln_VFT(T, A, B, Tinf):
      
      return np.log(A) + B / (T - Tinf)

# fractional linear solid?
