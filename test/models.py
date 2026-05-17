# updated version of this file is maintained at
# https://github.com/shullgroup/QKBPy/blob/main/test/models.py

import numpy as np
from scipy.optimize import curve_fit

# universal gas constant
R = 8.3145

# gaussian with baseline
def Gaussian(x, ctr, amp, wid, baseline=0):
        
        return baseline + amp * np.exp(-((x - ctr)/(2 * wid))**2)

def fitGaussian(df, xprop, yprop, **kwargs):
    '''
    Fits data to single Gaussian peak and adds to plot

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing DSC data.
    ax : mpl.axes.Axes
        Axes object to plot Gaussian fit on.
    return_err : bool, default False
        Option to return uncertainties for Tg and dT

    Returns
    -------
    ctr : float
        Center of Gaussian fit
    ctr_err : float
        Uncertainty in center of Gaussian fit
    wid : float
        Width of Gaussian fit
    wid_err : float
        Uncertainty in width of Gaussian fit

    '''
    
    ax = kwargs.get('ax', None)
    return_err = kwargs.get('return_err', False)
    guess = kwargs.get('guess', [100, 100, 100, 0])
    bounds = kwargs.get('bounds', ([0,1,1,-1000],[1e5,1e5,1e4,1000]))
    sigma = kwargs.get('sigma', None)
    absolute_sigma = kwargs.get('absolute_sigma', False)

    # clean up dataframe    
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df = df.dropna()
    
    #guesses for Temp, amplitude, and width for peak
    ctr_guess = df[yprop].max()
    guess = kwargs.get('guess', [df.query(f'{yprop} == @ctr_guess')[xprop].iloc[0], 1000, 1000, 0])
    
    #fit function to data and plot peak
    popt, pcov = curve_fit(Gaussian, df[xprop], df[yprop], p0=guess,
                           bounds=bounds, sigma=sigma, 
                           absolute_sigma=absolute_sigma, maxfev=5000)
    fit_x = np.linspace(df[xprop].min(),df[xprop].max(),num=1000)
    fit_y = Gaussian(fit_x, *popt)
    ctr = popt[0]
    wid = popt[2]
    
    # parameter uncertainties
    perr = np.sqrt(np.diag(pcov))
    ctr_err = perr[0]
    wid_err = perr[2]

    # add to a plot if axis given
    if ax:
        ax.plot(fit_x, fit_y, ':', color='k',
                label=f'Center = {ctr:0.0f} \n Width = {wid:0.0f}')
        ax.legend()
    
    return ctr, ctr_err, wid, wid_err

# arrhenius
def Arrhenius(T, A, Ea):

    return A * np.exp(Ea / (R * T))

# vft
def VFT(T, A, B, Tinf):
      
      return A * np.exp(B / (T - Tinf)) 

# ln versions
def ln_Arrhenius(T, A, Ea):
      
      return np.log(A) + Ea / (R * T)

def ln_VFT(T, A, B, Tinf):
      
      return np.log(A) + B / (T - Tinf)

def PowerLaw(x, A, n):

      return A * x**n

def ln_PowerLaw(x, A, n):
      
      return np.log(A) + n * x

# fractional linear solid?
