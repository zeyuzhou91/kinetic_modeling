"""
Functions related to arterial data and their processing. 
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from collections.abc import Callable
from numpy.typing import NDArray
from .tool import aux


# An array of numbers with time stamps
class TimedPoints:
    def __init__(self, 
                 name: str, 
                 data: NDArray, 
                 unit: str, 
                 t: NDArray, 
                 t_unit: str):
        
        self.name = name   # name of the quantity
        
        self.t = t     # time stamps
        self.data = data     # values of quantity
        
        self.t_unit = t_unit     # unit of time
        self.unit = unit    # unit of measured quantity
        
        self.f_to_fit = None   # Callable: time+params -> float
        self.f_fitted = None   # Callable: time -> float
        self.fitfuncname = None  # str, name of the fit function
        self.fit_params = None   # numpy.ndarray, parameters of the fitting function
        
    
    @classmethod
    def from_file(cls, filepath: str, name: str):
        """
        File must be a csv file with two columns. 
        """
        
        t, t_header, t_unit, data, header, unit = aux.read_from_csv_twocols(filepath)
        t = np.array(t)
        data = np.array(data)
        
        return cls(name, data, unit, t, t_unit)
        

    def plot(self, xlim: list[float] = None, ylim: list[float] = None) -> None:
        """
        Plot the data and fitted function. 
        """
        
        plt.figure()
        if self.t != [] and self.data != []: 
            plt.scatter(self.t, self.data, c='blue', label='data')
        if self.t != [] and self.f_fitted != None and self.fit_params.all() != None:
            tmax = np.max(self.t)
            tfit = np.linspace(0, tmax*1.1, 1000)
            yfit = self.f_fitted(tfit)
            plt.plot(tfit, yfit, c='red', label='fit')
        plt.xlabel(f't ({self.t_unit})')
        plt.ylabel(f'{self.y_unit}')
        plt.title(f'{self.name}')
        
        if xlim == None:
            pass
        else:
            plt.xlim(xlim)
        if ylim == None:
            pass
        else:
            plt.ylim(ylim)
            
        plt.legend()
        plt.show()
        
        return None
    
    
    def fit(self,
            f: Callable[..., float],
            fname: str,
            bounds: tuple | None = None, 
            ) -> None:
        """
        Fit the points to a parameterized function. 

        Parameters
        ----------
        f : function to be fitted, first argument of f must be time, the remaining arguments are parameters 
        fname : name of f 
        bounds : (optional) bounds on the parameters of f, see
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html for format
        

        """
        
        self.fitfuncname = fname
        self.f_to_fit = f
        
        if bounds is None:
            self.fit_params, _ = curve_fit(f, self.t, self.data)
        else:
            self.fit_params, _ = curve_fit(f, self.t, self.data, bounds=bounds)
            
        def f_fitted(t):
            return self.f_to_fit(t, *self.fit_params)
        self.f_fitted = f_fitted
        
        return None
    
    
    def print_fitparams(self):
        """
        Print the fitting parameters. 
        """
        
        print('Fitting parameters:\n')
        print(self.fitparams)
        return None




class BloodInput:
    def __init__(self,
                 ptac: TimedPoints,
                 pif: TimedPoints,
                 wb2p: TimedPoints,
                 frame_mid_points: NDArray,
                 unit: str,
                 t_unit: str):
        
        self.ptac = ptac
        self.pif = pif
        self.wb2p = wb2p
        
        self.frame_mid_points = frame_mid_points
        self.unit = unit
        self.t_unit = t_unit
        
        self.CP = lambda t: ptac.f_fitted(t) * pif.f_fitted(t)
        self.CB = lambda t: ptac.f_fitted(t) * wb2p.f_fitted(t)
        
        
    @classmethod
    def from_file(cls,
                  ptac_path: str,
                  pif_path: str,
                  p2wb_ratio_path: str,
                  frame_mid_points: NDArray,
                  ) -> None:
        
        # plasma concentration tac (before metabolite correction)
        ptac = TimedPoints.from_file(ptac_path, 'Plasma Activity Concentration')
        
        if ptac.unit == 'kBq/mL':
            pass
        elif ptac.unit == 'Bq/mL':
            ptac.data = ptac.data / 1000.0
            ptac.unit = 'kBq/mL'
        
        ptac.fit(f=zero_linear_3exp, fname='zero_linear_3exp')
        
        # plasma intact fraction
        pif = TimedPoints.from_file(pif_path, 'Plasma Intact Fraction')
        pif.fit(f=Hill, fname='Hill', bounds=Hill_bounds())
        
        
        # wb2plasma concentration ratio
        p2wb = TimedPoints.from_file(p2wb_ratio_path, 'Ratio of plasma-to-wholeblood activity concentration')
        
        wb2p = TimedPoints(name = 'Ratio of wholeblood-to-plasma activity concentration',
                           data = 1.0 / p2wb.data,
                           unit = p2wb.unit,
                           t = p2wb.t,
                           t_unit = p2wb.t_unit)
        wb2p.fit(f=twoExp, fname='twoExp', bounds=twoExp_bounds())
        
        
        return cls(ptac = ptac, 
                   pif = pif, 
                   wb2p = wb2p, 
                   frame_mid_points = frame_mid_points,
                   unit = ptac.unit,
                   t_unit = ptac.t_unit)
        
    

    def plot(self, trange: list[float]) -> None:
        """
        len(trange) = 2
        """

        ts = np.linspace(trange[0], trange[1], 1000)
        cps = self.CP(ts)
        cbs = self.CB(ts)
                
        plt.figure(1)
        plt.plot(ts, cps)
        plt.xlabel(f't ({self.t_unit})')
        plt.ylabel(f'{self.unit}')
        plt.title(r'$C_P(t)$)')
        plt.show()
        
        plt.figure(2)
        plt.plot(ts, cbs)
        plt.xlabel(f't ({self.t_unit})')
        plt.ylabel(f'{self.unit}')
        plt.title(r'$C_B(t)$)')
        plt.show()

        return None



# class BloodInput_old:
#     def __init__(self,
#                  CP: Callable[[float], float],
#                  CB: Callable[[float], float],
#                  frame_mid_points: NDArray,
#                  unit: str,
#                  t_unit: str):
        
#         self.CP = CP
#         self.CB = CB
#         self.frame_mid_points = frame_mid_points
#         self.unit = unit
#         self.t_unit = t_unit
        
        
#     def plot(self, trange: list[float]) -> None:
#         """
#         len(trange) = 2
#         """

#         ts = np.linspace(trange[0], trange[1], 1000)
#         cps = self.CP(ts)
#         cbs = self.CB(ts)
                
#         plt.figure(1)
#         plt.plot(ts, cps)
#         plt.xlabel(f't ({self.t_unit})')
#         plt.ylabel(f'{self.unit}')
#         plt.title(r'$C_P(t)$)')
#         plt.show()
        
        
#         plt.figure(2)
#         plt.plot(ts, cbs)
#         plt.xlabel(f't ({self.t_unit})')
#         plt.ylabel(f'{self.unit}')
#         plt.title(r'$C_B(t)$)')
#         plt.show()
        

#         return None




def zero_linear_3exp(t: float | np.ndarray, 
                     a, b, Tpk, A1, lamb1, A2, lamb2, lamb3,
                     ) -> float | np.ndarray:
    """
    Implements the following piecewise function f(t):
        
    f(t) = 0      if  t < -b/a
         = a*t+b  if -b/a <= t < Tpk
         = A1*exp{-lamb1*(t-Tpk)} + A2*exp{-lamb2*(t-Tpk)} + A3*exp{-lamb3*(t-Tpk)}  if t >= Tpk
    
    The parameters should satisfy:
        a * Tpk + b = A1 + A2 + A3
    
    Source: http://www.turkupetcentre.net/petanalysis/input_fitting_exp.html
    """
    
    A3 = a * Tpk + b - (A1 + A2)
    return (
        ((t>=0)&(t<-b/a)) * 0  +  
        ((t>=-b/a)&(t<Tpk)) * (a * t + b) + 
        (t>=Tpk) * (A1*np.exp(-lamb1*(t-Tpk)) + A2*np.exp(-lamb2*(t-Tpk)) + A3*np.exp(-lamb3*(t-Tpk))) )



def Hill(t: float | np.ndarray, a, b, c) -> float | np.ndarray:
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
    """
    
    return (1 - (1-a)*t**b/(c+t**b))


def Hill_bounds() -> tuple:
    """
    Returns bounds on the parameters a, b, c of the Hill function. 
    """
    
    # See details in the Hill function description
    bounds = ([0, 1,      0,    ],
              [1, np.inf, np.inf])
    
    return bounds
                


def twoExp(t, r1, r2, r3, r4) -> float:
    """
    Sum of two exponential functions. Often used for fitting plasma-to-blood 
    or blood-to-plasma ratio curves. 
    
    Constraints for all parameters: (0, +inf)
    
    Source: http://www.turkupetcentre.net/petanalysis/input_blood-to-plasma_fitting.html
    """
    
    return r1*np.exp(-r3*t) + r2*(1-np.exp(-r4*t))



def twoExp_bounds() -> tuple:
    """
    Returns bounds on the parameters r1, r2, r3, r4 of the twoExp function. 
    """
    
    bounds = ([0,      0,      0,    0],
              [np.inf, np.inf, np.inf, np.inf])
    
    return bounds


# oneExp is currently not used, maybe useful later
def oneExp(t, rmin, rmax, rate) -> float:
    """
    One exponential function f(t). 

    When t = 0, f(t) = rmax
    When t = +inf, f(t) = rmin    
    """

    return rmin + (rmax-rmin)*np.exp(-rate*t)


def oneExp_bounds() -> tuple:
    """
    Returns bounds on the parameters rmin, rmax, rate of the oneExp function. 
    """
    
    bounds = ([0,      0,      0],
              [np.inf, np.inf, np.inf])
    
    return bounds



# def read_plasma_tac(filepath: str) -> TimedPoints:
#     """
#     Read arterial plasma TAC information from a given file.
    
#     NOTE: this is before metabolite correction. 
#     """
    
#     # Initialize the plasma tac object
#     ptac = TimedPoints('Plasma Activity Concentration')
    
#     # Read the plasma tac information from file
#     t_data, t_header, t_unit, pconc_data, pconc_header, pconc_unit = aux.read_from_csv_twocols(filepath)
    
#     # Update plasma tac time attributes
#     ptac.t_data = t_data
#     ptac.t_unit = t_unit
    
    
#     if pconc_unit == 'kBq/mL':
#         ptac.y_data = pconc_data
#         ptac.y_unit = 'kBq/mL'
#     elif pconc_unit == 'Bq/mL':
#         # convert the pconc_data to unit [kBq/mL]
#         ptac.y_data = [x/1000.0 for x in pconc_data]
#         ptac.y_unit = 'kBq/mL'
    
#     return ptac
            


# def read_plasma_intact_frac(filepath: str) -> TimedPoints:
#     """
#     Read arterial plasma intact fraction information from a given file. 
#     """
    
#     # Initialize the object
#     pif = TimedPoints('Plasma Intact Fraction')
    
#     # Read the plasma intact fraction information from file
#     t_data, t_header, t_unit, pif_data, pif_header, pif_unit = aux.read_from_csv_twocols(filepath)
    
#     # Update time attributes
#     pif.t_data = t_data
#     pif.t_unit = t_unit
    
#     # Update fraction data
#     pif.y_data = pif_data
#     pif.y_unit = pif_unit  # should be "unitless"
    
#     return pif



# def read_wb2p_ratio(p2wb_filepath: str) -> TimedPoints:
#     """
#     Read wholeblood-to-plasma activity concentration ratio.  
    
#     NOTE: the input is p2wb, the result is the reverse: wb2p
#     """
    
#     # Initialize the object
#     wb2p = TimedPoints('Ratio of Wholeblood-to-plasma activity concentration')
    
#     # Read the data from file
#     t_data, t_header, t_unit, p2wb_data, p2wb_header, p2wb_unit = aux.read_from_csv_twocols(p2wb_filepath)
    
#     # Update time attributes
#     wb2p.t_data = t_data
#     wb2p.t_unit = t_unit
    
#     # Update fraction data
#     wb2p.y_data = [1.0/x for x in p2wb_data]
#     wb2p.y_unit = p2wb_unit     # should be unitless
    
#     return wb2p



# def generate_arterial_funcs(ptac_path: str,
#                             pif_path: str,
#                             p2wb_ratio_path: str,
#                             frameschedule: FrameSchedule,
#                             plot_details : bool 
#                             ) -> (BloodInput):
#     """
#     Generate two arterial functions, CP and CB. 
    
#     Parameters
#     ----------
#     ptac_path : file path of the plasma tac (before metabolite correction)
#     pif_path : file path of plasma intact/parent fraction
#     p2wb_ratio_path : file path of the plasma-to-wholeblood concentration ratio
#     plot_details : 
    
#     Outputs
#     ----------
#     CP : metabolitec-corrected plasma activity concentration curve, the arterial input function
#          unit: kBq/mL
#     CB : blood activity concentration curve
#          unit: kBq/mL
#     """
    
#     # plasma concentration tac (before metabolite correction)
#     ptac = read_plasma_tac(ptac_path)
#     ptac.fit(f=zero_linear_3exp, fname='zero_linear_3exp')
    

#     # plasma intact fraction
#     pif = read_plasma_intact_frac(pif_path)
#     pif.fit(f=Hill, fname='Hill', bounds=Hill_bounds())
    
    
#     # wb2plasma concentration ratio
#     wb2p = read_wb2p_ratio(p2wb_ratio_path)
#     wb2p.fit(f=twoExp, fname='twoExp', bounds=twoExp_bounds())
    
    
#     if plot_details:
#         ptac.plot(xlim=[0, 15])
#         pif.plot()
#         wb2p.plot()
    
    
#     binput = BloodInput(CP = lambda t: ptac.f_fitted(t) * pif.f_fitted(t), 
#                         CB = lambda t: ptac.f_fitted(t) * wb2p.f_fitted(t), 
#                         frameschedule = frameschedule,
#                         unit = ptac.y_unit,
#                         t_unit = ptac.t_unit)
    
    
#     return binput
    
    


if __name__ == "__main__":
    
    pass


        
        
    
    
    