"""
Auxiliary functions. 
"""

import os

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