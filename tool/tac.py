"""
Functions related to time-activity curves (TAC).
"""


import nibabel as nib
import numpy as np
import os
from . import aux  
from ..core import ROI

            


def extract_ROI_tac_from_PETimg(
        PETimg_path: str, 
        ROImask_path: str, 
        roi: ROI, 
        op_dir: str,
        PETimg_unit: str | None = None) -> str:
    """
    Extract ROI tac from the PET image and the corresponding mask. The roi object's 
    attributes will be updated. 
    
    Parameters
    ----------
    PETimg_path : file path. Path of the PET image, ending in .nii or .nii.gz
                  Can be either 3D (one frame) or 4D (multi-frame).
    ROImask_path : file path. Path of the ROI mask in PET domain, ending in .nii or .nii.gz
    roi : Contains information about the ROI. 
    op_dir : directory path. Path of the output directory, where the output file should be stored. 
    PETimg_unit : (optional) unit of the PET image. 

    Returns
    -------
    opfile_path : file path. Path of the output file. 

    """
    
    # For a given PET image, apply a binary mask and generate the ROI information. 
    
    # Load the input PET image
    PETimg = nib.load(PETimg_path)
    PETimg_data = PETimg.get_fdata()
    
    # Load the mask
    ROImask = nib.load(ROImask_path)
    ROImask_data = ROImask.get_fdata()
    
    # create the ROI object
    roi.num_voxels = np.count_nonzero(ROImask_data) 
    
    if len(PETimg.shape) == 3:
        # one frame (i.e. single 3D image)
        
        roi.num_frames = 1
        
        assert PETimg.shape == ROImask.shape, f"""The input image {PETimg_path} (shape = {PETimg_data.shape}) and 
        the mask {ROImask_path} (shape = mask_data.shape) should have the same dimension."""
        
        # outputs an array containing only the masked voxel values in PETimg
        ROI_reduced = np.extract(ROImask_data, PETimg_data)
        
        roi.tot_intensity = np.append(roi.tot_intensity, np.sum(ROI_reduced))
        roi.avg_intensity = np.append(roi.avg_intensity, np.mean(ROI_reduced))
        
        
    elif len(PETimg.shape) == 4:
        # multi-frame 
        
        num_frames = PETimg.shape[-1]
        roi.num_frames = num_frames
        
        assert PETimg.shape[0:3] == ROImask.shape, f"""Each frame of the input image {PETimg_path} (shape = {PETimg.shape[0:3]}) and 
        the mask {ROImask_path} (shape = ROImask_data.shape) should have the same dimension."""
    
        for i in range(num_frames):
            
            # the ith frame
            frame = PETimg_data[..., i]
            
            # outputs an array containing only the masked voxel values in this frame
            frame_ROI_reduced = np.extract(ROImask_data, frame)
            
            roi.tot_intensity = np.append(roi.tot_intensity, np.sum(frame_ROI_reduced))
            roi.avg_intensity = np.append(roi.avg_intensity, np.mean(frame_ROI_reduced))
            
           
    opfile_name = f'{roi.name}_avgIntensity.csv'
    opfile_path = os.path.join(op_dir, opfile_name)
    
    if PETimg_unit == 'Bq/mL':
        roi.avg_intensity /= 1000.0   # from Bq/mL to kBq/mL
            
    aux.write_to_csv_onecol(roi.avg_intensity, 'Avg Intensity', 'kBq/mL', opfile_path)
    
    return opfile_path
            
            


def extract_ROI_tac_from_csv(filepath: str, roi: ROI) -> None:
    """
    Extract ROI tac from a given csv file. The roi object's attributes will be updated. 
    
    Parameters
    ----------
    filepath : csv file path. 
    roi : the roi object, contains Contains information about the target ROI.  
    """
    
    tac, header, unit = aux.read_from_csv_onecol(filepath)
    roi.num_frames = len(tac)
    
    if unit == 'kBq/mL':
        roi.avg_intensity = np.array(tac)
    elif unit == 'Bq/mL':
        roi.avg_intensity = np.array(tac) / 1000.0
    
    return None
    

            