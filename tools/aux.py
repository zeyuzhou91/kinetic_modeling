"""
Auxiliary classes and functions.
"""

import os
import csv



class Environment:
    def __init__(self):
        
        self.root_dir = ''  # root directory path
        self.subj_dir = ''  # directory path of the subject
        self.MRI_dir = ''   # directory path for MRI images and segmentations
        self.PET_dir = ''   # directory path for PET images
        self.masks_dir = '' # directory path for masks 
        self.seg_path = ''  # path of the segmentation file
        self.mr2pet_lta_path = '' # path of the MR to PET linear transformation file
        self.pet2mr_lta_path = '' # path of the PET to MR linear transformation file
        self.framedurationfile_path = '' # path of the frame durations file


class ROI:
    def __init__(self, name):
        self.name = name     # string
        self.ID = None       # integer
        self.num_frames = None   # positive integer, 1 for single frame
        self.num_voxels = None  # integer
        self.vol_ml = None   # volume in [mL]
        
        self.avg_intensity = []  # list of decimals, average intensity of the voxels (mean of voxels) 
        self.tot_intensity = []  # list of decimals, total intensity of the voxels (sum of voxels)
        self.concentration = []     # list of decimals, average concentration in [Bq/mL], for dynamic imaging                           


class FrameSchedule:
    def __init__(self, durations):
        self.durations = durations   # list of numbers
        self.start_points = []       # list of numbers
        self.mid_points = []         # list of numbers
        self.calculate_attributes()
        
    def calculate_attributes(self):
        cur_start = 0.0
        for duration in self.durations:
            self.start_points.append(cur_start)
            
            nxt_start = cur_start + duration
            
            mid = (cur_start + nxt_start) / 2.0
            self.mid_points.append(mid)
            
            cur_start = nxt_start
    
        return None
            


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
    


def write_to_csv_onecol(header, arr, csvfile_path):
    """
    Write a 1D array to a CSV file as a column. 

    Parameters
    ----------
    header : string
    arr : list of numbers
    csvfile_path : string, file path, ending in .csv

    Returns
    -------
    None.

    """
    
    with open(csvfile_path, 'w', newline='') as csvfile:
        # Create a CSV writer object
        csv_writer = csv.writer(csvfile)

        # Write the header
        csv_writer.writerow([header])

        # Write the data
        for x in arr:
            csv_writer.writerow([x])
    
    return None

    

def write_to_csv_twocols(header1, arr1, header2, arr2, csvfile_path):
    """
    Write two 1D arrays to a CSV file as two columns.  

    Parameters
    ----------
    header1 : string
    arr1 : list of numbers
    header2 : string
    arr2 : list of numbers
    csvfile_path : string, file path, ending in .csv

    Returns
    -------
    None.

    """
    
    # Writing to CSV file
    with open(csvfile_path, 'w', newline='') as csvfile:
        # Create a CSV writer object
        csv_writer = csv.writer(csvfile)
    
        # Write the header
        csv_writer.writerow([header1, header2])
    
        # Write the data
        for x, y in zip(arr1, arr2):
            csv_writer.writerow([x, y])
    
    return None



def read_from_csv_onecol(csvfile_path):
    """
    Read a CSV file with one column, output the header and the data as a list. 

    Parameters
    ----------
    csvfile_path : string, file path, ending in .csv
        Only has one column, with the first row being the header. 

    Returns
    -------
    header : string
    data : list of floats
    """
    
    with open(csvfile_path, 'r') as csvfile:
        # Create a CSV reader object
        csv_reader = csv.reader(csvfile)

        # Read the header
        header = next(csv_reader)[0]

        # Read the data into a list
        data = [float(row[0]) for row in csv_reader]

    return header, data
    



            
if __name__ == "__main__":
    
    # test
    
    # file_path = "/home/documents/mask_cerebellum.nii.gz"
    
    # base, extension = extract_file_name(file_path)
    # print(base)
    # print(extension)
    
    # env = Environment()
    # env.MRI_dir = '/zeyu/documents/Data/FEPPA/F301/Freesurfer/bert1/mri'
    # env.PET_dir = '/zeyu/documents/Data/FEPPA/F301/PET'
    # print(env.MRI_dir)
    # print(env.PET_dir)
    
    header = 'concentration'
    arr = [10, 20, 30, 40, 50]
    write_to_csv_onecol(header, arr, '/Users/zeyuzhou/Documents/kinetic_modeling_test/FEPPA_20190523_AA_31002478/test.csv')
    
    
    # dur = [1, 2, 5, 6, 10]
    # FS = FrameSchedule(dur)
    # print(FS.durations)
    # print(FS.start_points)
    # print(FS.mid_points)
    
    pass
    
    
    
    