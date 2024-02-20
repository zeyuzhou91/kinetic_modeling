"""
Functions related to arterial data and their processing. 
"""

import numpy as np
from scipy.optimize import curve_fit
#import aux
from . import aux
#from aux import TimeCurve
#import aux.TimeCurve as TimeCurve




def zero_linear_3exp(t, a, b, Tpk, A1, lamb1, A2, lamb2, lamb3):
    """
    Implements the following function f(t):
        
    f(t) = 0      if  t < -b/a
         = a*t+b  if -b/a <= t < Tpk
         = A1*exp{-lamb1*(t-Tpk)} + A2*exp{-lamb2*(t-Tpk)} + A3*exp{-lamb3*(t-Tpk)}  if t >= Tpk
    
    The parameters should satisfy:
        a * Tpk + b = A1 + A2 + A3
    
    Source: http://www.turkupetcentre.net/petanalysis/input_fitting_exp.html

    Returns
    -------
    Float

    """
    
    #A3 = a * Tpk + b - (A1 + A2)
    return (
        ((t>=0)&(t<-b/a)) * 0  +  
        ((t>=-b/a)&(t<Tpk)) * (a * t + b) + 
        (t>=Tpk) * (A1*np.exp(-lamb1*(t-Tpk)) + A2*np.exp(-lamb2*(t-Tpk)) + (a*Tpk+b-(A1+A2))*np.exp(-lamb3*(t-Tpk))) )



def Hill(t, a, b, c):
    """
    Implements the Hill function f(t):
    
    f(t) = 1 - (1-a)t^b/(c+t^b)
    
    Constraints:
        0 <= a <= 1
        b >= 1
        c > 0
        
    This is a decreasing curve with values between 0 and 1. 
    Eventually, it stays at constant level a
    
    Source: http://www.turkupetcentre.net/petanalysis/input_parent_fitting_hill.html#:~:text=Hill%20type%20functions%20have%20been,of%20parent%20radiotracer%20in%20plasma.&text=%2C%20where%200%20%E2%89%A4%20a%20%E2%89%A4,parameter%20a%20the%20final%20level.
    
    There is a version of the Hill function with two more parameters d and e
    (see Source above), maybe useful for other situations. 
    
    
    Returns
    -------
    Float

    """
    
    return (1 - (1-a)*t**b/(c+t**b))




def fit_plasma_tac(tc):
    """
    Fit the arterial plasma TAC. 

    Parameters
    ----------
    tc : TimeCurve object
        Contains info regarding the arterial plasma TAC. 
    """
    
    tc.fitfunc = zero_linear_3exp
    tc.fitparams, _ = curve_fit(tc.fitfunc, tc.t_data, tc.y_data)
    
    return None
                


def fit_plasma_intact_frac(tc):
    """
    Fit the sampled arterial plasma intact fractions. 

    Parameters
    ----------
    tc : TimeCurve object
        Contains info regarding the sampled arterial plasma intact fractions.
    """
    
    tc.fitfunc = Hill
    
    # constraints of the Hill function parameters a, b, c
    # See details in the Hill function description
    bounds = ([0, 1,      0,    ],
              [1, np.inf, np.inf])
    
    tc.fitparams, _ = curve_fit(tc.fitfunc, tc.t_data, tc.y_data, bounds=bounds)
    
    return None


def fit_p2wb_ratio(p2wb):
    """
    Fit the sampled wholeblood-to-plasma concentration ratio.  

    Parameters
    ----------
    wb2p : TimeCurve object
        Contains info regarding the sampled wholeblood-to-plasma concentration ratio. 
    """
    
    
    # def f(t, r1, r2, r3, r4):
    #     return r1*np.exp(-r3*t) + r2*(1-np.exp(-r4*t))
        
    
    def f(t, rmin, rmax, rate):
        return rmin + (rmax-rmin)*np.exp(-rate*t)
    
    
    p2wb.fitfunc = f
    
    # bounds = ([0,      0,      0,    0],
    #           [np.inf, np.inf, np.inf, np.inf])
    
    bounds = ([0,      0,      0],
              [np.inf, np.inf, np.inf])
    
    p2wb.fitparams, _ = curve_fit(p2wb.fitfunc, p2wb.t_data, p2wb.y_data, bounds=bounds)
    
    return None





def read_plasma_tac(filepath):
    """
    Read arterial plasma TAC information and generate an object. 

    Parameters
    ----------
    filepath : string, file path
        Path of the arterial plasma TAC file. 

    Returns
    -------
    p_tac : TimeCurve object
        Contains info regarding the arterial plasma TAC. 

    """
    
    # Initialize the plasma tac object
    p_tac = aux.TimeCurve('plasma')
    
    # Read the plasma tac information from file
    t_header, t_unit, t_data, pconc_header, pconc_unit, pconc_data = aux.read_from_csv_twocols(filepath)
    
    # Update plasma tac time attributes
    p_tac.t_data = t_data
    p_tac.t_unit = t_unit
    
    if pconc_unit == 'Bq/mL':
        # convert the pconc_data to unit [kBq/mL]
        p_tac.y_data = [x/1000.0 for x in pconc_data]
        p_tac.y_unit = 'kBq/mL'
    
    return p_tac
            


def read_plasma_intact_frac(filepath):
    """
    Read arterial plasma intact fraction information and generate an object. 

    Parameters
    ----------
    filepath : string, file path
        Path of the arterial plasma intact fraction file. 

    Returns
    -------
    pif : TimeCurve object
        Contains info regarding the arterial plasma intact fraction.  

    """
    
    # Initialize the object
    pif = aux.TimeCurve('Plasma Intact Fraction')
    
    # Read the plasma intact fraction information from file
    t_header, t_unit, t_data, pif_header, pif_unit, pif_data = aux.read_from_csv_twocols(filepath)
    
    # Update plasma tac time attributes
    pif.t_data = t_data
    pif.t_unit = t_unit
    
    # Update fraction data
    pif.y_data = pif_data
    pif.y_unit = pif_unit
    
    return pif



def read_plasma2wb_ratio(filepath):
    """
    Read wholeblood-to-plasma activity concentration ratio. 
    
    Parameters
    ----------
    filepath : string, file path
        Path of the file containing the plasma-to-wholeblood ratio. 
        NOTE: This is NOT wholeblood-to-plasma

    Returns
    -------
    wb2p : TimeCurve object
        Contains info regarding the wholeblood-to-plasma activity concentration ratio.  
    """
    
    # Initialize the object
    p2wb = aux.TimeCurve('Plasma-to-wholeblood activity concentration ratio')
    
    # Read the plasma intact fraction information from file
    t_header, t_unit, t_data, p2wb_header, p2wb_unit, p2wb_data = aux.read_from_csv_twocols(filepath)
    
    # Update plasma tac time attributes
    p2wb.t_data = t_data
    p2wb.t_unit = t_unit
    
    # Update fraction data
    p2wb.y_data = p2wb_data
    p2wb.y_unit = p2wb_unit
    
    return p2wb



def generate_AIF(ptac_path, pif_path, p2wb_conc_ratio_path):
    
    # plasma concentration tac
    ptac = read_plasma_tac(ptac_path)
    fit_plasma_tac(ptac)
    #ptac.plot(xlim=[0, 15])
    #ptac.print_fitparams()
    
    # plasma intact fraction
    pif = read_plasma_intact_frac(pif_path)
    fit_plasma_intact_frac(pif)
    #pif.plot()
    #pif.print_fitparams()
    
    
    # plasms2wb concentration ratio
    p2wb = read_plasma2wb_ratio(p2wb_conc_ratio_path)
    fit_p2wb_ratio(p2wb)
    #p2wb.plot(ylim=[1.0, 1.8])
    #p2wb.print_fitparams()
    
    def AIF(t):
        
        return ptac.fitfunc(t, *ptac.fitparams) * pif.fitfunc(t, *pif.fitparams)
    
    
    def WB_tac(t):
        
        # TO UPDATE
        return ptac.fitfunc(t, *ptac.fitparams)
    
    
    return AIF, WB_tac
    
    


if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
    

    
    ptac_path = '/Users/zeyuzhou/Documents/kinetic_modeling/test/plasma_tac.csv'
    pif_path = '/Users/zeyuzhou/Documents/kinetic_modeling/test/plasma_intact_fraction.csv'
    p2wb_conc_ratio_path = '/Users/zeyuzhou/Documents/kinetic_modeling/test/plasma2wb_conc_ratio_VAT025.csv'
    
    AIF, WB_tac = generate_AIF(ptac_path, pif_path, p2wb_conc_ratio_path)
    
    ts = np.linspace(0, 160, 1000)
    ys = AIF(ts)
    plt.figure()
    plt.plot(ts, ys, c='red')
    plt.xlabel('Time (min)')
    plt.ylabel('kBq/mL')
    plt.xlim([0, 15])
    #plt.legend()
    plt.show()
        
        
    
    
    