"""
Auxiliary classes and functions.
"""

import os



class Environment:
    def __init__(self):
        
        self.MRI_dir = ''   # directory path for MRI images and segmentations
        self.PET_dir = ''   # directory path for PET images
        self.masks_dir = '' # directory path for masks 
        self.seg_path = ''  # path of the segmentation file
        self.mr2pet_lta_path = '' # path of the MR to PET linear transformation file
        self.pet2mr_lta_path = '' # path of the PET to MR linear transformation file



class ROI:
    def __init__(self, name):
        self.name = name     # string
        self.ID = None       # integer
        self.numvoxels = None  # integer
        self.vol_ml = None   # volume in [mL]
        self.avg_intensity = None  # average intensity of the voxels (mean of voxels) 
        self.tot_intensity = None  # total intensity of the voxels (sum of voxels)
        self.conc = None    # decimal, average concentration in [Bq/mL]
        self.num_frames = None   # positive integer, 1 for single frame
        self.concs = []     # list of decimals, average concentration in [Bq/mL], for dynamic imaging
                            



def extract_file_name(file_path):
    """
    Extract the file base name and extension from its path. 
    
    Example: /home/documents/mask_cerebellum.nii.gz
    Return: mask_cerebellum and .nii.gz

    Parameters
    ----------
    file_path : string, file path
        The full path of the file.

    Returns
    -------
    basename : string
        The base name of the file. 
    extension: string
        The extension of the file. 
    """

    # os.path.basename returns the full name of the file, including the extension
    fullname = os.path.basename(file_path)
    
    # split fullname from the first occurrence of '.'
    # e.g.: mask_cerebellum.nii.gz => mask_cerebellum and nii.gz
    basename, extension = fullname.split('.', 1)
    extension = '.' + extension
    
    return basename, extension
    
    
    

            
if __name__ == "__main__":
    
    # test
    file_path = "/home/documents/mask_cerebellum.nii.gz"
    
    base, extension = extract_file_name(file_path)
    print(base)
    print(extension)
    
    env = Environment()
    env.MRI_dir = '/zeyu/documents/Data/FEPPA/F301/Freesurfer/bert1/mri'
    env.PET_dir = '/zeyu/documents/Data/FEPPA/F301/PET'
    print(env.MRI_dir)
    print(env.PET_dir)
    