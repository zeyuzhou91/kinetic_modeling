"""
Functions related to time-activity curves (TAC).
"""

import os
import nibabel as nib
import numpy as np
from numpy.typing import NDArray
from . import aux  


            


def extract_tac_from_PETimg(
        PETimg_path: str, 
        ROImask_path: str,
        ROIname: str,
        op_dir: str,
        opfile_name: str | None = None,
        PETimg_unit: str | None = None) -> (NDArray, str):
    """
    Extract ROI tac from the PET image and the corresponding mask. The roi object's 
    attributes will be updated. 
    
    Parameters
    ----------
    PETimg_path : file path. Path of the PET image, ending in .nii or .nii.gz
                  Can be either 3D (one frame) or 4D (multi-frame).
    ROImask_path : file path. Path of the ROI mask in PET domain, ending in .nii or .nii.gz
    ROIname : name of the ROI
    op_dir : directory path. Path of the output directory, where the output file should be stored. 
    opfile_name : (optional) name of the output filename. 
    PETimg_unit : (optional) unit of the PET image. 

    Returns
    -------
    avg_intensity : average intensity of the ROI at each frame
    unit : output unit of the tac

    """
    
    # For a given PET image, apply a binary mask and generate the ROI information. 
    
    # Load the input PET image
    PETimg = nib.load(PETimg_path)
    PETimg_data = PETimg.get_fdata()
    
    # Load the mask
    ROImask = nib.load(ROImask_path)
    ROImask_data = ROImask.get_fdata()
    
    # num_voxels = np.count_nonzero(ROImask_data) 
    
    if len(PETimg.shape) == 3:
        # one frame (i.e. single 3D image)
        
        assert PETimg.shape == ROImask.shape, f"""The input image {PETimg_path} (shape = {PETimg_data.shape}) and 
        the mask {ROImask_path} (shape = mask_data.shape) should have the same dimension."""
        
        # outputs an array containing only the masked voxel values in PETimg
        ROI_reduced = np.extract(ROImask_data, PETimg_data)
                
        tot_intensity = np.array([np.sum(ROI_reduced)])
        avg_intensity = np.array([np.mean(ROI_reduced)])
        
        
    elif len(PETimg.shape) == 4:
        # multi-frame 
        
        num_frames = PETimg.shape[-1]
        
        assert PETimg.shape[0:3] == ROImask.shape, f"""Each frame of the input image {PETimg_path} (shape = {PETimg.shape[0:3]}) and 
        the mask {ROImask_path} (shape = ROImask_data.shape) should have the same dimension."""
    
        tot_intensity = np.zeros(num_frames)
        avg_intensity = np.zeros(num_frames)        
        for i in range(num_frames):
            
            # the ith frame
            frame = PETimg_data[..., i]
            
            # outputs an array containing only the masked voxel values in this frame
            frame_ROI_reduced = np.extract(ROImask_data, frame)
            
            tot_intensity[i] = np.sum(frame_ROI_reduced)
            avg_intensity[i] = np.mean(frame_ROI_reduced)
            
           
    if opfile_name is None: 
        opfile_name = f'{ROIname}_avgIntensity.csv'
    opfile_path = os.path.join(op_dir, opfile_name)
    
    if PETimg_unit == 'kBq/mL':
        aux.write_to_csv_onecol(avg_intensity, 'Avg Intensity', 'kBq/mL', opfile_path)
        unit = PETimg_unit
    elif PETimg_unit == 'Bq/mL':
        avg_intensity /= 1000.0   # from Bq/mL to kBq/mL
        aux.write_to_csv_onecol(avg_intensity, 'Avg Intensity', 'kBq/mL', opfile_path)
        unit = 'kBq/mL'
    elif PETimg_unit == 'unitless':
        aux.write_to_csv_onecol(avg_intensity, 'Avg Intensity', 'unitless', opfile_path)
        unit = 'unitless'
    
    return avg_intensity, unit
            
            


def extract_tac_from_csv(filepath: str) -> (NDArray, str):
    """
    Extract ROI tac from a given csv file. The tac0's attributes will be updated. 
    
    Parameters
    ----------
    filepath : csv file path. 
    
    Returns
    -------
    ys : tac values 
    unit : unit of ys
    
    """
    
    ys, header, unit = aux.read_from_csv_onecol(filepath)
    ys = np.array(ys)
    
    if unit == 'kBq/mL':
        pass

    elif unit == 'Bq/mL':
        ys = ys / 1000.0
        unit = 'kBq/mL'

    elif unit == 'unitless':
        pass
        
    return ys, unit
    

            