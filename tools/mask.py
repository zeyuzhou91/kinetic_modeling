"""
Functions related to masking. 

"""

import numpy as np
import nibabel as nib
import copy
import os



def foo():
    print("Hello, this is the foo function in mask.py.")
    return None


def create_MR_mask_including(in_IDs, segfile_path, opROI_name, op_dir):
    """
    Create an MR mask (binary .nii.gz image) that includes the given ROIs. 

    Parameters
    ----------
    in_IDs : a list of integers
        The integer IDs of ROIs that the mask should include. 
    segfile_path : string 
        The path of the MR segmentation file from which the mask is created. 
        The file should be .nii or .nii.gz 
    opROI_name : string 
        The name of the output combined ROI. 
    op_dir : string
        The path of the output directory where the output mask is stored. 

    Returns
    -------
    opfile_name: string
        The name of the output mask file, ending in .nii.gz
    opfile_path: string
        The path of the output mask file. 

    """
    
    # load segmentation as an nifti image
    seg = nib.load(segfile_path)
    
    # make a copy of the image (_data is np array)
    mask_opROI_data = copy.deepcopy(seg.get_fdata())
    
    # Find the binary mask for the output ROI, including the in_IDs
    mask_opROI_data = np.isin(mask_opROI_data, in_IDs).astype(float)
    
    # Make the nifti mask image
    mask_opROI = nib.Nifti1Image(mask_opROI_data, seg.affine, seg.header)
    
    opfile_name = f'mask_mr_{opROI_name}.nii.gz'
    opfile_path = os.path.join(op_dir, opfile_name)
    
    # Save the nifti mask image
    nib.save(mask_opROI, opfile_path)
    
    return opfile_name, opfile_path
    
    


def create_MR_mask_excluding(ex_IDs, segfile_path, opROI_name, op_dir):
    """
    Create an MR mask (binary .nii.gz image) that excludes the given ROIs. 

    Parameters
    ----------
    ex_IDs : a list of integers
        The integer IDs of ROIs that the mask should exclude. 
    segfile_path : string 
        The path of the MR segmentation file from which the mask is created. 
        The file should be .nii or .nii.gz 
    opROI_name : string 
        The name of the output combined ROI. 
    op_dir : string
        The path of the output directory where the output mask is stored. 

    Returns
    -------
    opfile_name: string
        The name of the output mask file, ending in .nii.gz
    opfile_path: string
        The path of the output mask file. 

    """
    
    # load segmentation as an nifti image
    seg = nib.load(segfile_path)
    
    # make a copy of the image (_data is np array)
    mask_opROI_data = copy.deepcopy(seg.get_fdata())
    
    # Find the binary mask for the output ROI, excluding the ex_IDs
    mask_opROI_data = (~np.isin(mask_opROI_data, ex_IDs)).astype(float)
    
    # Make the nifti mask image
    mask_opROI = nib.Nifti1Image(mask_opROI_data, seg.affine, seg.header)
    
    opfile_name = f'mask_mr_{opROI_name}.nii.gz'
    opfile_path = os.path.join(op_dir, opfile_name)
    
    # Save the nifti mask image
    nib.save(mask_opROI, opfile_path)
    
    return opfile_name, opfile_path



            
            

            

            
            
            
            

            

            
        
            
            
            
            