"""
Core class definitions. 
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from .tool import aux


class Environment:
    def __init__(self):
        
        self.root_dir = ''  # root directory path
        self.subj_dir = ''  # directory path of the subject
        self.MRI_dir = ''   # directory path for MRI images and segmentations
        self.PET_dir = ''   # directory path for PET images
        self.tac_dir = ''   # directory path for ROI tac information
        self.AIF_dir = ''   # directory path for AIF data
        self.km_dir = ''    # directory path for kinetic modeling (temp)
        self.masks_dir = '' # directory path for masks 
        self.seg_path = ''  # path of the segmentation file
        self.mr2pet_lta_path = '' # path of the MR to PET linear transformation file
        self.pet2mr_lta_path = '' # path of the PET to MR linear transformation file
        self.framedurationfile_path = '' # path of the frame durations file
        self.framemidpointfile_path = '' # path of the frame mid-points file 
        self.AIF_ptac_path = ''
        self.AIF_pif_path = ''
        self.AIF_p2wb_ratio_path = ''


class ROI:
    def __init__(self, name):
        self.name = name     # str
        self.ID = None       # int
        self.num_frames = None   # positive int, 1 for single frame
        self.num_voxels = None  # int
        self.vol_ml = None   # volume in [mL]
        
        self.avg_intensity = np.array([])  # average intensity of the voxels (mean of voxels) 
        self.tot_intensity = np.array([])  # total intensity of the voxels (sum of voxels)
        self.concentration = np.array([])  # average concentration in [Bq/mL], for dynamic imaging          

         
        self.twotcm_params = {'K1': None, 'k2': None, 'k3': None, 'k4':None, 'VND':None, 'VT':None, 'VS':None, 'BPND':None}              



class FrameSchedule:
    def __init__(self, mid_points: NDArray, durations: NDArray, start_points: NDArray):
        self.mid_points = mid_points
        self.durations = durations         
        self.start_points = start_points             
        
    @classmethod
    def from_midpoints(cls, filepath: str):
        mid_points, _, __ = aux.read_from_csv_onecol(filepath)
        durations = []
        start_points = []
        
        cur_start = 0.0
        for mid in mid_points:
            
            start_points.append(cur_start)
            
            duration = 2 * (mid - cur_start)
            durations.append(duration)
            
            nxt_start = cur_start + duration
            
            cur_start = nxt_start
        
        mid_points = np.array(mid_points)
        durations = np.array(durations)
        start_points = np.array(start_points)
        
        return cls(mid_points, durations, start_points)

    
    @classmethod
    def from_durations(cls, filepath: str):
        durations, _, __ = aux.read_from_csv_onecol(filepath)
        mid_points = []
        start_points = []
        
        cur_start = 0.0
        for duration in durations:
            
            start_points.append(cur_start)
            
            nxt_start = cur_start + duration
            
            mid = (cur_start + nxt_start) / 2.0
            mid_points.append(mid)
            
            cur_start = nxt_start
        
        mid_points = np.array(mid_points)
        durations = np.array(durations)
        start_points = np.array(start_points)
        
        return cls(mid_points, durations, start_points)


            

# Function of time
class FuncOfTime:
    def __init__(self):
        self.f = None          # function, input: number, output: number
        self.t_unit = None     # str, unit of time
        self.y_unit = None     # str, unit of output
        self.name = ''         # str, name of function 
        
    def plot(self, trange: list[float]) -> None:
        """
        len(trange) = 2
        
        """
        
        assert self.f != None, "Function not defined yet."
        
        ts = np.linspace(trange[0], trange[1], 1000)
        ys = self.f(ts)
                
        plt.figure()
        plt.plot(ts, ys)
        plt.xlabel(f't ({self.t_unit})')
        plt.ylabel(f'{self.y_unit}')
        plt.title(f'{self.name}')
        plt.show()

        return None
    


    

        
            
if __name__ == "__main__":
        
    
    pass
    
    
    
    