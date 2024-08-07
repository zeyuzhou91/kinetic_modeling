"""
Functions related to masking. 

"""

import numpy as np
import nibabel as nib
import copy
import os
from skimage.morphology import ball, binary_dilation
from . import aux  # NOTE: simply "import aux" won't work


def create_MR_mask_including(
        in_IDs: list[int], 
        seg_path: str, 
        opROI_name: str, 
        op_dir: str) -> str:
    """
    Create an MR mask (binary .nii.gz image) that includes the given ROIs. 

    Parameters
    ----------
    in_IDs : The integer IDs of ROIs that the mask should include. 
    seg_path : file path. The path of the MR segmentation file from which the 
               mask is created. The file should be .nii or .nii.gz.
    opROI_name : The name of the output combined ROI. 
    op_dir : directory path. The path of the output directory where the output mask is stored. 

    Returns
    -------
    opmask_path: file path. The path of the output mask file, ending in .nii.gz.

    """
    
    # load segmentation as an nifti image
    seg = nib.load(seg_path)
    
    # make a copy of the image (_data is np array)
    mask_opROI_data = copy.deepcopy(seg.get_fdata())
    
    # Find the binary mask for the output ROI, including the in_IDs
    mask_opROI_data = np.isin(mask_opROI_data, in_IDs).astype(int)
    
    # Make the nifti mask image
    mask_opROI = nib.Nifti1Image(mask_opROI_data, seg.affine, seg.header)
    
    opmask_fullname = f'mask_mr_{opROI_name}.nii.gz'
    opmask_path = os.path.join(op_dir, opmask_fullname)
    
    # Save the nifti mask image
    nib.save(mask_opROI, opmask_path)
    
    return opmask_path
    
    


def create_MR_mask_excluding(
        ex_IDs: list[int], 
        seg_path: str, 
        opROI_name: str, 
        op_dir: str) -> str:
    """
    Create an MR mask (binary .nii.gz image) that excludes the given ROIs. 

    Parameters
    ----------
    ex_IDs : The integer IDs of ROIs that the mask should exclude. 
    seg_path : file path. The path of the MR segmentation file from which the 
               mask is created. The file should be .nii or .nii.gz.
    opROI_name : The name of the output combined ROI. 
    op_dir : directory path. The path of the output directory where the output mask is stored. 

    Returns
    -------
    opmask_path: file path. The path of the output mask file, ending in .nii.gz

    """
    
    # load segmentation as an nifti image
    seg = nib.load(seg_path)
    
    # make a copy of the image (_data is np array)
    mask_opROI_data = copy.deepcopy(seg.get_fdata())
    
    # Find the binary mask for the output ROI, excluding the ex_IDs
    mask_opROI_data = (~np.isin(mask_opROI_data, ex_IDs)).astype(int)
    
    # Make the nifti mask image
    mask_opROI = nib.Nifti1Image(mask_opROI_data, seg.affine, seg.header)
    
    opmask_fullname = f'mask_mr_{opROI_name}.nii.gz'
    opmask_path = os.path.join(op_dir, opmask_fullname)
    
    # Save the nifti mask image
    nib.save(mask_opROI, opmask_path)
    
    return opmask_path





def linear_transform(
        ipmask_path: str, 
        inDomain: str, 
        outDomain: str, 
        lta_path: str, 
        thr: float, 
        save_bfthr_mask: bool, 
        op_dir: str) -> str:
    """
    Performs linear transform of the input mask from inDomain to outDomain. 

    Parameters
    ----------
    ipmask_path : directory path. The path of the input mask .nii.gz file. 
    inDomain : The input domain, 'mr' or 'pet'.
    outDomain : The output domain, 'mr' or 'pet'.
    lta_path : file path. The path of the .reg.lta file, containing information of the linear transform. 
    thr : float in [0, 1]. The threshold for mapping decimal values to binary 
          values for the PET mask transformed from MR domain.
    save_bfthr_mask : True - save the intermediate decimal-valued mask before thresholding; 
                      False - do not save.  
    op_dir : directory path. The path of the output directory where the output mask is stored. 

    Returns
    -------
    opmask_path: file path. The path of the output mask file, ending in .nii.gz

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
    
    # Check if the lta_path exists
    assert os.path.isfile(lta_path) is True, f"""The lta_path {lta_path} does not exist."""
    
    # map the input mask from inDomain to outDomain
    # mapped mask has decimal values
    os.system(f'$FREESURFER_HOME/bin/mri_convert -at {lta_path} {ipmask_path} {opmask_bfthr_path}')
    
    # load the bfthr output mask
    opmask_bfthr = nib.load(opmask_bfthr_path)
    opmask_bfthr_data = opmask_bfthr.get_fdata()
    
    # thresholding the bfthr output mask (decimal) to make it binary
    opmask_data = (opmask_bfthr_data >= thr).astype(int)
    opmask = nib.Nifti1Image(opmask_data, opmask_bfthr.affine, opmask_bfthr.header)
    
    # name and path of the binary output mask
    opmask_fullname = opmask_basename + extension
    opmask_path = os.path.join(op_dir, opmask_fullname)
    nib.save(opmask, opmask_path)
        
    if save_bfthr_mask == False:
        # delete the intermediate bfthr mask
        os.remove(opmask_bfthr_path)

    return opmask_path
            



def create_PET_mask_including(
        in_IDs: list[int], 
        thr: float, 
        save_PET_bfthr_mask: bool, 
        save_MR_mask: bool, 
        opROI_name: str, 
        seg_path: str,
        mr_masks_dir: str,
        mr2pet_lta_path: str,
        op_dir: str) -> str:
    """
    Create a PET mask (binary .nii.gz image) that includes the given ROIs. 
    
    Parameters
    ----------
    in_IDs : The integer IDs of ROIs that the mask should include. 
    thr : float in [0, 1]. The threshold for mapping decimal values to binary 
          values for the PET mask transformed from MR domain.
    save_PET_bfthr_mask : True - save the intermediate decimal-valued PET mask before thresholding; 
                          False - do not save. 
    save_MR_mask : True - save the intermediate MR mask; 
                   False - do not save.
    opROI_name : The name of the output combined ROI. 
    seg_path: segmentation file path
    mr_masks_dir: directory path of the MR masks folder
    mr2pet_lta_path: file path of the mr2pet lta transformation
    op_dir : directory path. The path of the output directory where the output mask is stored.  

    Returns
    -------
    PETmask_path : file path of the PET mask.
    """

    MRmask_path = create_MR_mask_including(
        in_IDs = in_IDs,
        seg_path = seg_path,
        opROI_name = opROI_name,
        op_dir = mr_masks_dir)

    PETmask_path = linear_transform(
        ipmask_path = MRmask_path,
        inDomain = 'mr',
        outDomain = 'pet',
        lta_path = mr2pet_lta_path,
        thr = thr,
        save_bfthr_mask = save_PET_bfthr_mask,
        op_dir = op_dir)

    if save_MR_mask == False:
        # delete the MR mask
        os.remove(MRmask_path)
        
    return PETmask_path

       

def create_PET_mask_excluding(
        ex_IDs: list[int], 
        thr: float, 
        save_PET_bfthr_mask: bool, 
        save_MR_mask: bool, 
        opROI_name: str, 
        seg_path: str,
        mr_masks_dir: str,
        mr2pet_lta_path: str,
        op_dir: str) -> str:
    """
    Create a PET mask (binary .nii.gz image) that excludes the given ROIs. 
    
    Parameters
    ----------
    ex_IDs : The integer IDs of ROIs that the mask should exclude. 
    thr : float in [0, 1]. The threshold for mapping decimal values to binary 
          values for the PET mask transformed from MR domain.
    save_PET_bfthr_mask : True - save the intermediate decimal-valued PET mask before thresholding; 
                          False - do not save. 
    save_MR_mask : True - save the intermediate MR mask; 
                   False - do not save.
    opROI_name : The name of the output combined ROI. 
    seg_path: segmentation file path
    mr_masks_dir: directory path of the MR masks folder
    mr2pet_lta_path: file path of the mr2pet lta transformation
    op_dir : directory path. The path of the output directory where the output mask is stored. 

    Returns
    -------
    PETmask_path : file path of the PET mask.
    """

    MRmask_path = create_MR_mask_excluding(
        ex_IDs = ex_IDs,
        seg_path = seg_path,
        opROI_name = opROI_name,
        op_dir = mr_masks_dir)
        
    PETmask_path = linear_transform(
        ipmask_path = MRmask_path,
        inDomain = 'mr',
        outDomain = 'pet',
        lta_path = mr2pet_lta_path,
        thr = thr,
        save_bfthr_mask = save_PET_bfthr_mask,
        op_dir = op_dir)

    if save_MR_mask == False:
        # delete the MR mask
        os.remove(MRmask_path)
        
    return PETmask_path        




def generate_masked_img(ipimg_path, mask_path, maskedROI_name, op_dir):
    """
    For a given input image, apply a binary mask and generate a masked image. 

    Parameters
    ----------
    ipimg_path : string, file path
        The file path of the input image, ending in .nii or .nii.gz. 
        The input image can be either 3D (one frame) or 4D (multi-frame).
    mask_path : string, file path
        The file path of the binary mask, ending in .niior .nii.gz.
        The mask should be 3D. 
    maskedROI_name : string
        The name of the masked ROI.
    op_dir : string, directory path
        The path of the output directory where the masked image is stored. 

    Returns
    -------
    opimage_path : string, file path
        The path of the output masked image file, ending in .nii or .nii.gz.
        The output image is of the same dimension as the input image. 
    """
    
    # Load the input image
    ipimg = nib.load(ipimg_path)
    ipimg_data = ipimg.get_fdata()
    
    # Load the mask
    mask = nib.load(mask_path)
    mask_data = mask.get_fdata()
    

    
    if len(ipimg.shape) == 3:
        # one frame (i.e. single 3D image)
        
        assert ipimg.shape == mask.shape, f"""The input image {ipimg_path} (shape = {ipimg_data.shape}) and 
        the mask {mask_path} (shape = mask_data.shape) should have the same dimension."""
    
        # Generate the output nifti image
        opimg_data = ipimg_data * mask_data
        
    elif len(ipimg.shape) == 4:
        # multi-frame 
        
        num_frames = ipimg.shape[-1]
        
        assert ipimg.shape[0:3] == mask.shape, f"""Each frame of the input image {ipimg_path} (shape = {ipimg.shape[0:3]}) and 
        the mask {mask_path} (shape = mask_data.shape) should have the same dimension."""
    
        opimg_data = copy.deepcopy(ipimg_data)
        for i in range(num_frames):
            opimg_data[..., i] = ipimg_data[..., i] * mask_data


    opimg = nib.Nifti1Image(opimg_data, ipimg.affine, ipimg.header)


    ipimg_basename, extension = aux.extract_file_name(ipimg_path)
        
    # Generate the output image's name and path
    # E.g. if input image = frame5.nii.gz and maskedROI_name is "cerebellum"
    # Then output image = frame5_cerebellum.nii.gz
    opimg_basename = ipimg_basename + '_' + maskedROI_name
    opimg_fullname = opimg_basename + extension
    opimg_path = os.path.join(op_dir, opimg_fullname)


    nib.save(opimg, opimg_path)
    
    return opimg_path
            
        


def dilation(
        ball_radius: float,
        ippath: str,
        oppath: str) -> None:
    """
    Dilate the mask by a ball with a given radius.

    """

    # load input image
    ip = nib.load(ippath)
    
    # make a copy of the image 
    op_data = copy.deepcopy(ip.get_fdata())
    
    # dilation by a 
    op_data = binary_dilation(op_data, ball(ball_radius)).astype(op_data.dtype)
    
    # Make the output image
    op = nib.Nifti1Image(op_data, ip.affine, ip.header)
            
    # Save the output image
    nib.save(op, oppath)
    
    return None 




def union(
        ippath1: str,
        ippath2: str,
        oppath: str) -> None:
    """
    Find the union of two masks.
    """

    ip1 = nib.load(ippath1)
    ip1_data = copy.deepcopy(ip1.get_fdata()).astype(int)
    
    ip2 = nib.load(ippath2)
    ip2_data = copy.deepcopy(ip2.get_fdata()).astype(int)
    
    op_data = ip1_data | ip2_data
    
    # Make the output image
    op = nib.Nifti1Image(op_data, ip1.affine, ip1.header)
            
    # Save the output image
    nib.save(op, oppath)
    
    return None


def intersect(
        ippath1: str,
        ippath2: str,
        oppath: str) -> None:
    """
    Find the intersection of two masks.
    """

    ip1 = nib.load(ippath1)
    ip1_data = copy.deepcopy(ip1.get_fdata()).astype(int)
    
    ip2 = nib.load(ippath2)
    ip2_data = copy.deepcopy(ip2.get_fdata()).astype(int)
    
    op_data = ip1_data & ip2_data
    
    # Make the output image
    op = nib.Nifti1Image(op_data, ip1.affine, ip1.header)
            
    # Save the output image
    nib.save(op, oppath)
    
    return None


def minus(
        ippath1: str,
        ippath2: str,
        oppath: str) -> None:
    """
    Find all regions in image1 that is not in image2. 
    """

    ip1 = nib.load(ippath1)
    ip1_data = copy.deepcopy(ip1.get_fdata()).astype(int)
    
    ip2 = nib.load(ippath2)
    ip2_data = copy.deepcopy(ip2.get_fdata()).astype(int)
    
    intersect = ip1_data & ip2_data 
    
    op_data = np.logical_xor(ip1_data, intersect)

    # Make the output image
    op = nib.Nifti1Image(op_data, ip1.affine, ip1.header)
            
    # Save the output image
    nib.save(op, oppath)
    
    return None
