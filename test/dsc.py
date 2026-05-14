#Import any individual functions from outside packages 
#that are used in your functions.
#These are called dependencies.
# updated version of this file is maintained at
# https://github.com/shullgroup/rheoQCM/blob/master/QCMFuncs/DSC_functions.py
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy.signal import savgol_filter

from scipy.optimize import curve_fit
from .utils import *

def readDSC(path, **kwargs):
    '''
    Read txt file from DSC experiment and convert to a DataFrame
    
    Parameters
    ----------
    path : Path
        Path object to the .txt file containing the DSC data

    mode : str, default 'conv'
        Flag for mode of DSC, i.e. conventional (conv) vs. moduluated (mdsc)

    apply_savgol : bool, default True
        Option to apply the Savitsky-Golay filter, especially useful
        for the derivatives.

    savgol_window : int, default 151
        Window length for filtering derivative data.  Must be odd int

    savgol_polyoder : int, default 4
        Polynomial degree for filtering derivative data.  Must be 
        less than window length.
    
    Returns
    -------
    df : pd.DataFrame
        DataFrame containing relevant experimental data, i.e.
        temperature, norm. heat flow, deriv. heat flow
        
    '''
    mode = kwargs.get('mode', 'conv')
    apply_savgol = kwargs.get('apply_savgol', True)
    savgol_window = kwargs.get('savgol_window', 151)
    savgol_polyorder = kwargs.get('savgol_polyorder', 4)

    if mode == 'mdsc':
        start_row = first_line(path)
        with open(path, 'r') as f:
            df = pd.read_csv(f, delimiter="\t", skiprows=start_row,
                            usecols=[0,1,2,3,7],
                            names=['time','temp','q_rev','q_non','dq_revdT'])

        df = df.dropna()
        if apply_savgol:
            # smooth q_rev data
            # uses Savitsky-Golay filter with default window 151 and polyorder 4
            # because raw exported data extremely coarse
            df['q_rev'] = savgol_filter(df['q_rev'], 
                                    window_length=savgol_window, 
                                    polyorder=savgol_polyorder)
            df['dq_revdT'] = savgol_filter(df['dq_revdT'],
                                           window_length=savgol_window,
                                           polyorder=savgol_polyorder)


        return df
    
    elif mode == 'conv':
        start_row = first_line(path)
        with open(path, 'r') as f:
            df = pd.read_csv(f, delimiter="\t", skiprows=start_row,
                            usecols=[0,1,2],
                            names=['time','temp','q'])

        df = df.dropna()
        if apply_savgol:
            # calculate derivative of q w.r.t. T
            # uses Savitsky-Golay filter with default window 151 and polyorder 4
            # because derivative data extremely noisy/jumpy
            df['dqdT'] = savgol_filter(np.gradient(df['q'], df['temp']), 
                                    window_length=savgol_window, 
                                    polyorder=savgol_polyorder)
        else:
            df['dqdT'] = np.gradient(df['q'], df['temp'])


        return df

def readMDSC(path, **kwargs):
    """
    Reads a Modulated Differential Scanning Calorimetry (MDSC) data file and 
    processes it into a pandas DataFrame.
    
    Parameters:
    -----------
    path : str
        The file path to the MDSC data file in tab-separated value (TSV) format.
    
    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the processed MDSC data with the following 
        columns:
        - 't': Time converted to seconds (originally in minutes in the file).
        - 'temp': Temperature.
        - 'q_norm': Normalized heat flow.
        - 'q_rev': Reversible heat flow.
        - 'q_nrev': Non-reversible heat flow.
        - 'q_tot': Total heat flow.
        - 'cp_rev': Reversible heat capacity.
        - 'cpdT': Heat capacity derivative.
        
    Notes:
    ------
    - The input file should have at least 9 header rows which are skipped 
      during reading.
    - The data is expected to have 8 columns which are read and named 
      accordingly.
    - Rows with any NaN values are removed from the resulting DataFrame.
    - Time values are converted from minutes to seconds assuming the DSC 
      software uses minutes as the standard unit.
    """
    start_row = first_line(path)
    df = pd.read_csv(path, sep='\t', skiprows=start_row,
                     usecols =[0,1,2,3,4,5,6,7],
                     names = ['t', 'temp', 'q_norm', 'q_rev', 'q_nrev', 
                              'q_tot','cp_rev','dcpdT'])
    # remove rows with any nan values
    df = df.dropna()
    
    # convert time to seconds (assuming the dsc software always has minues
    # s the standard unit)
    df['t'] = 60 * df['t']
         
    return df


    
def plotDSC(df, **kwargs):
    '''
    Generate typical plots for DSC experiments. Emphasis primarily placed
    on finding Tg as opposed to other transitions for now.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing experimental data read in from the readDSC function
    mode : str, default 'conv'
        DSC mode used for experiment. Options are 'conv' for conventional DSC
        (i.e. heating and cooling ramps) or 'mdsc' for temperature modulated DSC.
    ax : mpl.axes.Axes, default None
        Axes for the heat flow if one already exists.
    twin : mpl.axes.Axes, default None
        Twin axis for the derivative data if one already exists.
    title : str, default None
        Title for the plot.
    savepath : Path, default None
        Path for saving the plot if you want to.
    min_temp : float, default None
        Minimum temperature to use for the DSC experiment.
    max_temp : float, default None
        Maximum temperature to use for the DSC experiment.
    showTg : bool, default True
        Option to show or not show Tg chosen from maximum in derivative
        e.g. if using the fitGaussianDSC function below to find Tg.
    no_legend : bool, default False
        Option to remove legend e.g. if using legend generated from fitGaussianDSC
    legendloc : int or str, default 0
        Location of legend if need to manually change
    legendsize : float, default 10.
        Size of legend if need to manually change
    show_deriv_ticks : bool, default False
        Option to show the tick marks on the twin axis for the derivative heat
        flow. Default False because typically too crowded.
    orientation : str, default 'exo_up'
        Orientation/convention for heat flow annotation. Typically 'Exo Up' but also
        accepts 'endo_up' for 'Endo Up' annotation.

    Returns
    -------
    Tg : float
        Glass transition temperature (Tg) in deg. Celsius chosen as maximum in 
        derivative of heat flow.
    ax : mpl.axes.Axes
        Axes instance used for plotting the heat flow.  Output used for carrying to
        other plots for overlays and such.
    twin : mpl.axes.Axes
        Axes instance used for plotting the derivative heat flow.  Output used
        for carrying to fitGaussianDSC function for better calculation of Tg

    '''

    mode = kwargs.get('mode', 'conv')
    ax = kwargs.get('ax', None)
    twin = kwargs.get('twin', None)
    title = kwargs.get('title', None)
    savepath = kwargs.get('savepath', None)
    min_temp = kwargs.get('min_temp', None)
    max_temp = kwargs.get('max_temp', None)
    no_legend = kwargs.get('no_legend', False)
    legendloc = kwargs.get('legendloc', 0)
    legendsize = kwargs.get('legendsize', 10.)
    show_deriv_ticks = kwargs.get('show_deriv_ticks', False)
    orientation = kwargs.get('orientation', 'exo_up')

    # replace infs with nan and drop nans
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df = df.dropna()
    
    # filter temperatures
    if min_temp != None:
        df = df.query('temp > @min_temp')
    else:
        min_temp = df['temp'].min()
    if max_temp != None:
        df = df.query('temp < @max_temp')
    else:
        max_temp = df['temp'].max()
    

    # choose which mode: conventional or mdsc
    if mode == 'conv':
        q = 'q'
        peak = df['dqdT'].min()
        Tg = float(df.query('dqdT == @peak')['temp'].iloc[0])
    elif mode == 'mdsc':
        q = 'q_rev'
        peak = df['dq_revdT'].min()
        Tg = float(df.query('dq_revdT == @peak')['temp'].iloc[0])

    # create ax for typical DSC plot with derivative to show Tg
    if not ax:
        fig, ax = plt.subplots(1,1, figsize=(4,3), constrained_layout=True)
        ax.set_xlabel('Temperature ($^\\circ$C)')
        ax.set_xlim([min_temp-5, max_temp+15])
        ax.set_ylabel('Norm. Heat Flow (W/g)')
        ax.set_title(title)
    # create twin for heat flow derivative
    if not twin:
        twin = ax.twinx()
        twin.set_ylabel('Deriv. Heat Flow (W/g*$^{\\circ}$C)')
    twin._get_lines = ax._get_lines

    # plot q and dqdT
    ax.plot(df['temp'], df[f'{q}'], '-', label='q')
    twin.plot(df['temp'], -df[f'd{q}dT'], '--', label='dq/dT')
    
    # fit data to Gaussian for better Tg selection
    Tg, dT = fitGaussianDSC(df, ax=twin)
    
    # make sure the legends all go together
    if no_legend:
        pass
    else:
        handles, labels = ax.get_legend_handles_labels()
        handles_t, labels_t = twin.get_legend_handles_labels()
        ax.legend(
            handles=[*handles,*handles_t],
            labels = [*labels,*labels_t],
            loc=legendloc,
            prop={'size':legendsize})
        
    # annotation for heat flow convention
    if orientation == 'exo_up':
        ax.annotate('Exo Up', (5,5), xycoords='axes points')
    elif orientation == 'endo_up':
        ax.annotate('Endo Up', (5,5), xycoords='axes points')

    # option to show the tick marks for the derivative
    # but usually they crowd the figure too much
    if not show_deriv_ticks:   
        twin.set_yticks([])

    # option to save the figure
    if savepath:
        plt.savefig(savepath)

    return Tg, ax, twin

def fitGaussianDSC(df, ax, **kwargs):
    '''
    Fits data to single Gaussian peak and adds to DSC plot

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing DSC data.
    ax : mpl.axes.Axes
        Axes object to plot Gaussian fit on. Typically twin from DSC plot.
    bounds : tuple, default ([-100,0,0,-1],[200,1,30,1])
        Bounds for the fitting parameters.
    guess : list, default [df.query('dqdT == @maxdqdT')['temp'].iloc[0],
                                 3e-2, 10, 0]
        Initial guess (p0) for fitting.

    Returns
    -------
    Tg : float
        Glass transition temperature from center of Gaussian fit in deg. C
    Tg_err : float
        Uncertainty in Tg fit in deg. C
    dT : float
        Breadth of glass transition from width of Gaussian fit in deg. C
    dT_err : float
        Uncertainty in dT fit in deg. C

    '''

    bounds = kwargs.get('bounds', ([-100,0,0,-1],[200,1,30,1]))

    def Gaussian(x, *params):
        y = np.zeros_like(x)
        ctr = params[0]
        amp = params[1]
        wid = params[2]
        y0 = params[3]
        y = y + y0 + amp*np.exp(-((x - ctr)/(2*wid))**2)
        return y

    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df = df.dropna()
    
    #guesses for Temp, amplitude, and width for peak
    maxdqdT = df['dqdT'].min()
    guess = kwargs.get('guess', [df.query('dqdT == @maxdqdT')['temp'].iloc[0],
                                 3e-2, 10, 0])
    
    #fit function to data and plot peak
    popt, pcov = curve_fit(Gaussian, df['temp'], -df['dqdT'], 
                           p0=guess, 
                           bounds=bounds)
    fit_temp = np.linspace(df['temp'].min(),df['temp'].max(),num=1000)
    fit = Gaussian(fit_temp, *popt)
    dT = popt[2]
    Tg = popt[0]
    # parameter uncertainties
    perr = np.sqrt(np.diag(pcov))
    Tg_err = perr[0]
    dT_err = perr[2]

    # plot the fit
    ax.plot(fit_temp, fit, ':', color='k',
            label=f'T$_g = {Tg:0.1f} ^\\circ$C \n $\u03b4T = {dT:0.1f} ^\\circ$C')
    
    return Tg, Tg_err, dT, dT_err
