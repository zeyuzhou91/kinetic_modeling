"""
Math functions. 

"""

import os
import numpy as np
import nibabel as nib
import copy



def multiply_local(
        value: int | float, 
        ipimg: str,
        opimg: str) -> None:
    """
    Multiply the image by a given value.
    
    Do it locally, read the local ipimg and produce the opimg to the same folder. 
    """
    
    # load input image
    ip = nib.load(ipimg)
    
    # make a copy of the image 
    op_data = copy.deepcopy(ip.get_fdata())
    
    # Multiply
    op_data = op_data * value
    
    # Make the output image
    op = nib.Nifti1Image(op_data, ip.affine, ip.header)
            
    # Save the output image
    nib.save(op, opimg)
    
    return None



def multiply(
        value: int | float, 
        ippath: str,
        oppath: str) -> None:
    """
    Multiply the image by a given value.
    """
    
    # load input image
    ip = nib.load(ippath)
    
    # make a copy of the image 
    op_data = copy.deepcopy(ip.get_fdata())
    
    # Multiply
    op_data = op_data * value
    
    # Make the output image
    op = nib.Nifti1Image(op_data, ip.affine, ip.header)
            
    # Save the output image
    nib.save(op, oppath)
    
    return None



def thresholding_local(
        lb: float,
        ub: float,
        ipimg: str,
        opimg: str) -> None:
    """
    Threshold the image by the given lower-bound and higher-bound. 

    """

    # load input image
    ip = nib.load(ipimg)
    
    # make a copy of the image 
    op_data = copy.deepcopy(ip.get_fdata())
    
    # thresholding
    op_data = np.clip(op_data, lb, ub)
    
    # Make the output image
    op = nib.Nifti1Image(op_data, ip.affine, ip.header)
            
    # Save the output image
    nib.save(op, opimg)
    
    return None
            
        
def thresholding(
        lb: float,
        ub: float,
        ippath: str,
        oppath: str) -> None:
    """
    Threshold the image by the given lower-bound and higher-bound. 

    """

    # load input image
    ip = nib.load(ippath)
    
    # make a copy of the image 
    op_data = copy.deepcopy(ip.get_fdata())
    
    # thresholding
    op_data = np.clip(op_data, lb, ub)
    
    # Make the output image
    op = nib.Nifti1Image(op_data, ip.affine, ip.header)
            
    # Save the output image
    nib.save(op, oppath)
    
    return None   
            

def max_value_in_image(ippath: str) -> float:
    """
    Find the max value in a given image. 

    """

    # load input image
    ip = nib.load(ippath)
    
    op_data = ip.get_fdata()
    
    max_v = np.max(op_data)
    
    return max_v


def percentile_value_in_image(q: float, 
                              ippath: str) -> float:
    """
    Find the q-th percentile of all pixel/voxel values in an image. 
    
    q: in [0, 100]

    """

    # load input image
    ip = nib.load(ippath)
    
    op_data = ip.get_fdata()
    
    r = np.percentile(op_data.flatten(), q)
    
    return r



        
def add(
        ippath1: str,
        ippath2: str,
        oppath: str) -> None:
    """
    Add two images. They must be of the same dimension. 

    """
    
    ip1 = nib.load(ippath1)
    ip1_data = copy.deepcopy(ip1.get_fdata())
    
    ip2 = nib.load(ippath2)
    ip2_data = copy.deepcopy(ip2.get_fdata())
    
    op_data = ip1_data + ip2_data
    
    # Make the output image
    op = nib.Nifti1Image(op_data, ip1.affine, ip1.header)
            
    # Save the output image
    nib.save(op, oppath)
    
    return None
            