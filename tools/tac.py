"""
Functions related to time-activity curves (TAC).
"""


import nibabel as nib
import numpy as np
import os
from . import aux  # NOTE: simply "import aux" won't work

            

# Function: 
# Input: a multi-frame PET image, a target ROI (name or index), the frame schedule
# Output: a csv file with two columes, one for frame time (in mins), the other for ROI activity concentration (Bq/mL)
# Use functions in mask.py



def extract_PET_ROI_info(PETimg_path, ROImask_path, ROI, op_dir):
    """

    Parameters
    ----------
    PETimg_path : string, file path
        Path of the PET image, ending in .nii or .nii.gz
        Can be either 3D (one frame) or 4D (multi-frame).
    ROImask_path : string, file path
        Path of the ROI mask in PET domain, ending in .nii or .nii.gz
    ROI : ROI object
        Contains information about the ROI. 
    op_dir : string, directory path
        Path of the output directory, where the output file should be stored. 

    Returns
    -------
    opfile_path : string, file path
        Path of the output file. 

    """
    
    # For a given PET image, apply a binary mask and generate the ROI information. 
    
    # Load the input PET image
    PETimg = nib.load(PETimg_path)
    PETimg_data = PETimg.get_fdata()
    
    # Load the mask
    ROImask = nib.load(ROImask_path)
    ROImask_data = ROImask.get_fdata()
    
    # create the ROI object
    ROI.num_voxels = np.count_nonzero(ROImask_data) 
    
    if len(PETimg.shape) == 3:
        # one frame (i.e. single 3D image)
        
        ROI.num_frames = 1
        
        assert PETimg.shape == ROImask.shape, f"""The input image {PETimg_path} (shape = {PETimg_data.shape}) and 
        the mask {ROImask_path} (shape = mask_data.shape) should have the same dimension."""
        
        # outputs an array containing only the masked voxel values in PETimg
        ROI_reduced = np.extract(ROImask_data, PETimg_data)
        
        ROI.tot_intensity.append(np.sum(ROI_reduced))
        ROI.avg_intensity.append(np.mean(ROI_reduced))
        
        
    elif len(PETimg.shape) == 4:
        # multi-frame 
        
        num_frames = PETimg.shape[-1]
        ROI.num_frames = num_frames
        
        assert PETimg.shape[0:3] == ROImask.shape, f"""Each frame of the input image {PETimg_path} (shape = {PETimg.shape[0:3]}) and 
        the mask {ROImask_path} (shape = ROImask_data.shape) should have the same dimension."""
    
        for i in range(num_frames):
            
            # the ith frame
            frame = PETimg_data[..., i]
            
            # outputs an array containing only the masked voxel values in this frame
            frame_ROI_reduced = np.extract(ROImask_data, frame)
            
            ROI.tot_intensity.append(np.sum(frame_ROI_reduced))
            ROI.avg_intensity.append(np.mean(frame_ROI_reduced))
            
           
    opfile_name = f'{ROI.name}_avgIntensity.csv'
    opfile_path = os.path.join(op_dir, opfile_name)
    aux.write_to_csv_onecol('Avg Intensity', ROI.avg_intensity, opfile_path)
    
    return opfile_path
            
            




            