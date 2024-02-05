"""
Functions related to masking. 

"""

import numpy as np
import nibabel as nib
import copy
import os
from . import aux  # NOTE: simply "import aux" won't work



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
    segfile_path : string, file path
        The path of the MR segmentation file from which the mask is created. 
        The file should be .nii or .nii.gz 
    opROI_name : string 
        The name of the output combined ROI. 
    op_dir : string, directory path
        The path of the output directory where the output mask is stored. 

    Returns
    -------
    opmask_path: string, file path
        The path of the output mask file, ending in .nii.gz

    """
    
    # load segmentation as an nifti image
    seg = nib.load(segfile_path)
    
    # make a copy of the image (_data is np array)
    mask_opROI_data = copy.deepcopy(seg.get_fdata())
    
    # Find the binary mask for the output ROI, including the in_IDs
    mask_opROI_data = np.isin(mask_opROI_data, in_IDs).astype(float)
    
    # Make the nifti mask image
    mask_opROI = nib.Nifti1Image(mask_opROI_data, seg.affine, seg.header)
    
    opmask_fullname = f'mask_mr_{opROI_name}.nii.gz'
    opmask_path = os.path.join(op_dir, opmask_fullname)
    
    # Save the nifti mask image
    nib.save(mask_opROI, opmask_path)
    
    return opmask_path
    
    


def create_MR_mask_excluding(ex_IDs, segfile_path, opROI_name, op_dir):
    """
    Create an MR mask (binary .nii.gz image) that excludes the given ROIs. 

    Parameters
    ----------
    ex_IDs : a list of integers
        The integer IDs of ROIs that the mask should exclude. 
    segfile_path : string, file path 
        The path of the MR segmentation file from which the mask is created. 
        The file should be .nii or .nii.gz 
    opROI_name : string 
        The name of the output combined ROI. 
    op_dir : string, directory path
        The path of the output directory where the output mask is stored. 

    Returns
    -------
    opmask_path: string, file path
        The path of the output mask file, ending in .nii.gz

    """
    
    # load segmentation as an nifti image
    seg = nib.load(segfile_path)
    
    # make a copy of the image (_data is np array)
    mask_opROI_data = copy.deepcopy(seg.get_fdata())
    
    # Find the binary mask for the output ROI, excluding the ex_IDs
    mask_opROI_data = (~np.isin(mask_opROI_data, ex_IDs)).astype(float)
    
    # Make the nifti mask image
    mask_opROI = nib.Nifti1Image(mask_opROI_data, seg.affine, seg.header)
    
    opmask_fullname = f'mask_mr_{opROI_name}.nii.gz'
    opmask_path = os.path.join(op_dir, opmask_fullname)
    
    # Save the nifti mask image
    nib.save(mask_opROI, opmask_path)
    
    return opmask_path




            
            
def linear_transform(ipmask_path, inDomain, outDomain, lta_path, thr, save_bfthr_mask, op_dir):
    """
    Performs linear transform of the input mask from inDomain to outDomain. 

    Parameters
    ----------
    ipmask_path : string, directory path
        The path of the input mask .nii.gz file. 
    inDomain : string
        The input domain, 'mr' or 'pet'
    outDomain : string
        The output domain, 'mr' or 'pet'
    lta_path : string, file path
        The path of the .reg.lta file, containing information of the linear transform. 
    thr : float in [0, 1]
        The threshold for deciding 0 or 1 for decimal voxel values in output mask. 
    save_bfthr_mask : boolean
        True: save the intermediate decimal-valued mask before thresholding; False: do not save.  
    op_dir : string, directory path
        The path of the output directory where the output mask is stored. 

    Returns
    -------
    opmask_path: string, file path
        The path of the output mask file, ending in .nii.gz

    """
    
    # dashed versions of inDomain and outDomain
    inD = f'_{inDomain}_'
    outD = f'_{outDomain}_'
    
    ipmask_basename, extension = aux.extract_file_name(ipmask_path)

    assert inD in ipmask_basename, f"The mask name {ipmask_basename} should contain {inD}"
    
    # create the output mask'a base name by replacing the first occurence of inD in ipmask_name with outD
    opmask_basename = ipmask_basename.replace(inD, outD, 1)
    
    # names and path of the bfthr output mask
    opmask_bfthr_basename = opmask_basename + '_bfthr'
    opmask_bfthr_fullname = opmask_bfthr_basename + extension
    opmask_bfthr_path = os.path.join(op_dir, opmask_bfthr_fullname)
    
    # map the input mask from inDomain to outDomain
    # mapped mask has decimal values
    os.system(f'$FREESURFER_HOME/bin/mri_convert -at {lta_path} {ipmask_path} {opmask_bfthr_path}')
    
    # load the bfthr output mask
    opmask_bfthr = nib.load(opmask_bfthr_path)
    opmask_bfthr_data = opmask_bfthr.get_fdata()
    
    # thresholding the bfthr output mask (decimal) to make it binary
    opmask_data = (opmask_bfthr_data >= thr).astype(float)
    opmask = nib.Nifti1Image(opmask_data, opmask_bfthr.affine, opmask_bfthr.header)
    
    # name and path of the binary output mask
    opmask_fullname = opmask_basename + extension
    opmask_path = os.path.join(op_dir, opmask_fullname)
    nib.save(opmask, opmask_path)
        
    if save_bfthr_mask == False:
        # delete the intermediate bfthr mask
        os.remove(opmask_bfthr_path)

    return opmask_path
            
            

            
#def generate_masked_img(ipimage_path, mask_path)
            
        
            
            
            
            