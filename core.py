"""
Core class definitions. 
"""

import os
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from collections.abc import Callable
from .tool import aux, mask, tac



class Environment:
    def __init__(self, path_names_file):
        
        self.path_dict = {}
        
        self.read_path_names_from_file(path_names_file)
        
        
    def print_all(self):
        
        for (name, path) in self.path_dict.items():
            print(f'{name} = {path}')
        
        return None
    
    
    def reset_all(self):
        
        for name in self.path_dict.keys():
            self.path_dict[name] = ''
            
        
    def read_path_names_from_file(self, file_path):
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    # Strip whitespace and add non-empty lines to the array
                    path_name = line.strip()
                    if path_name:
                        self.path_dict[path_name] = ''
        except FileNotFoundError:
            print(f"Error: File not found at {file_path}")
        except IOError:
            print(f"Error: Could not read file at {file_path}")
        
        return None



class ROI:
    def __init__(self, name, IDs=None, ID_type=None):
        self.name = name     # str
        self.IDs = IDs if IDs is not None else None     # list[int]
        self.ID_type = ID_type if ID_type is not None else None  # "including" or "excluding"

        # self.num_voxels = None  # int
        # self.vol_ml = None   # volume in [mL]
                 

                    


class FrameSchedule:
    def __init__(self, mid_points: NDArray, durations: NDArray, start_points: NDArray, unit: str):
        self.mid_points = mid_points
        self.durations = durations         
        self.start_points = start_points
        self.unit = unit
        self.num_frames = len(mid_points)
        
    @classmethod
    def from_midpoints(cls, 
                       array: list[float] = None, 
                       filepath: str = None):
    
        if array is not None and filepath is not None:
            raise ValueError("Only one of array or filepath should be provided, not both.")
        
        elif array is not None:
            mid_points = array
            unit = ''
            
        elif filepath is not None:
            mid_points, _, unit = aux.read_from_csv_onecol(filepath)
        
        else:
            raise ValueError("Either array or filepath must be provided.")
        
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
        
        return cls(mid_points, durations, start_points, unit)

    
    @classmethod
    def from_durations(cls, 
                       array: list[float] = None,
                       filepath: str = None):
        
        if array is not None and filepath is not None:
            raise ValueError("Only one of array or filepath should be provided, not both.")
        
        elif array is not None:
            durations = array
            unit = ''
            
        elif filepath is not None:
            durations, _, unit = aux.read_from_csv_onecol(filepath)
        
        else:
            raise ValueError("Either array or filepath must be provided.")
        
        
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
        
        return cls(mid_points, durations, start_points, unit)


            

class TAC:
    def __init__(self, 
                 frameschedule: FrameSchedule, 
                 data: NDArray | None = None,
                 rois: list[ROI] | None = None,
                 unit: str | None = None):
        
        self.frameschedule = frameschedule
        self.data = data
        self.rois = rois
        self.unit = unit
        self.num_elements = None
        
        if data is not None:
            if data.shape[1] != self.frameschedule.num_frames:
                raise ValueError("data.shape[1] must be the same as num_frames in frameschedule")
            
        if data is not None and rois is not None:
            if data.shape[0] != len(rois):
                raise ValueError("data.shape[0] must equal to len(rois)")
            self.num_elements = data.shape[0]
            
    
    def read_rois_from_file(self, file_path):
        self.rois = []
        
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    if line.strip():
                        parts = line.strip().split(' ')
                        name = parts[0]
                        ID_type = parts[1]
                        IDs_str = parts[2:]
                        IDs = [int(x) for x in IDs_str]
                        roi = ROI(name = name,
                                  IDs = IDs,
                                  ID_type = ID_type)
                        self.rois.append(roi)
                        
                        self.num_elements = len(self.rois)
        except FileNotFoundError:
            print(f"Error: File not found at {file_path}")
        except IOError:
            print(f"Error: Could not read file at {file_path}")
            
        return None
    

    def extract_tac(self, 
                    tacfile_name_extension: str,
                    PETimg_path: str | None = None, 
                    PETimg_unit: str | None = None,
                    env: Environment | None = None) -> None:
        
        """
        Extract tac from the tac file (csv) if it exists; otherwise, extract it 
        from the PET image. 
        
        PETimg_path: dynamic (4D) image. 
        """
        
        if self.rois is None:
            raise ValueError("self.rois not defined")
        
        self.data = np.zeros((self.num_elements, self.frameschedule.num_frames))
        
        
        for (i, roi) in enumerate(self.rois):
        
            tacfile_name = roi.name + tacfile_name_extension
            tacfile_path = os.path.join(env.path_dict['tacs_dir'], tacfile_name)
            
            if os.path.exists(tacfile_path):
                
                self.data[i,:], self.unit = tac.extract_tac_from_csv(tacfile_path)
                
            else:
                
                THR = 0.8  # adjustable
                            
                if roi.ID_type == "including":
                
                    ROImask_path = mask.create_PET_mask_including(
                        in_IDs = roi.IDs,
                        thr = THR,
                        save_PET_bfthr_mask = False,
                        save_MR_mask = True,
                        opROI_name = roi.name, 
                        seg_path = env.path_dict['seg_path'],
                        mr_masks_dir = env.path_dict['mr_masks_dir'],
                        mr2pet_lta_path = env.path_dict['mr2pet_lta_path'],
                        op_dir = env.path_dict['pet_masks_dir'])
                    
                elif roi.ID_type == "excluding":
                    
                    ROImask_path = mask.create_PET_mask_excluding(
                        in_IDs = roi.IDs,
                        thr = THR,
                        save_PET_bfthr_mask = False,
                        save_MR_mask = True,
                        opROI_name = roi.name, 
                        seg_path = env.path_dict['seg_path'],
                        mr_masks_dir = env.path_dict['mr_masks_dir'],
                        mr2pet_lta_path = env.path_dict['mr2pet_lta_path'],
                        op_dir = env.path_dict['pet_masks_dir'])
    
                self.data[i,:], self.unit = tac.extract_tac_from_PETimg(
                    PETimg_path = PETimg_path,
                    ROImask_path = ROImask_path,
                    ROIname = roi.name,
                    op_dir = env.path_dict['tacs_dir'],
                    PETimg_unit = PETimg_unit)
        
        
        if self.data.shape[1] != self.frameschedule.num_frames:
            raise ValueError("self.data.shape[1] must be the same as num_frames in frameschedule")
        
        
        return None


       
    # def plot(self) -> None:

    #     # check if self.frameschedule and self.ys exist before plotting

    #     ts = self.frameschedule.mid_points
    #     ys = self.ys
        
    #     plt.figure()
    #     plt.scatter(ts, ys)
    #     plt.xlabel(f't ({self.frameschedule.unit})')
    #     plt.ylabel(f'{self.unit}')
    #     if self.ROI is not None:
    #         plt.title(f'{self.ROI.name}')
    #     plt.show()

    #     return None
    
    




# # Function of time
# class FuncOfTime:
#     def __init__(self):
#         self.f = None          # function, input: number, output: number
#         self.t_unit = None     # str, unit of time
#         self.y_unit = None     # str, unit of output
#         self.name = ''         # str, name of function 
        
#     def plot(self, trange: list[float]) -> None:
#         """
#         len(trange) = 2
        
#         """
        
#         assert self.f != None, "Function not defined yet."
        
#         ts = np.linspace(trange[0], trange[1], 1000)
#         ys = self.f(ts)
                
#         plt.figure()
#         plt.plot(ts, ys)
#         plt.xlabel(f't ({self.t_unit})')
#         plt.ylabel(f'{self.y_unit}')
#         plt.title(f'{self.name}')
#         plt.show()

#         return None
    


    

        
            
if __name__ == "__main__":
        
    
    array = [1, 1, 1, 2, 2, 3, 3, 4, 6, 10]
    myfs = FrameSchedule.from_durations(array)
    myfs.unit = 'min'
    
    mytac = TAC(myfs)
    
    print(mytac.frameschedule.durations)

    
    # myf = lambda t: t**2
    
    # ctac = ContinuousTAC(frameschedule = myfs, 
    #                      f = myf, 
    #                      unit = 'Bq/mL', 
    #                      name = 'ctac')
    
    # ctac.plot(trange = [0, 2])
    
    
    
    ys = [2, 3, 5, 6, 8, 9, 10, 13, 15, 19]
    roi = ROI(name = 'brain', IDs = [1, 2], ID_type = "including")
    
    dtac = DiscreteTAC(frameschedule = myfs,
                       ys = ys, 
                       unit = 'Bq/mL', 
                       roi = roi)
    
    dtac.plot()
    
    
    
    
    
    