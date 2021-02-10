#  Uses  python/3.6.4-anaconda.
#  Originally written by Bin Yang, and adapted by Kevin Dean.
#  UTSW microscope software saves data in the form of Cell1/1_CH(channel)_(time).tif
#  Code is designed to parse each Cell#, and identify the number of channels (channel) and timepoints.
#  Works on Linux, would need to be modified for Windows OS.
#  Written by a biochemist, so use at your own risk and expect bugs.
#
#  python deskewDirectory.py /path/to/your/directory.

import os
import numpy as np # 1.18.1
import sys
import re
from tifffile import imread, imsave
from multiprocessing import Pool

#  Specify the number of parallel processes you would like to use.  Depends on your filesize and available RAM.
numberThreads = 15

# Specify the size of your Z-step (the scan direction).
dz = 800

# Specify your lateral pixel size.
xypixelsize = 115

# Specify the angle of your oblique illumination.
angle = 30

#  Parse the Command Line Inputs
parent_directory = sys.argv[1]

#  Confirm that the string concludes with a forward slash
if parent_directory.endswith('/'):
    print(parent_directory)
else:
    print('parent directory needs /')
    parent_directory = parent_directory + '/'
    print(parent_directory)

def parse_directory(parent_directory):
    #  Determine number of Cell# subdirectories in the parent directory
    path_contents = os.listdir(parent_directory)
    loop_idx = 0
    cell_number = 0

    #  Iterate through each experiment.  Cell 1, Cell 2, Cell3...
    for cellIdx in path_contents:
        if "Cell" in path_contents[cell_number]:

            #  Iterate through the contents of each experiment
            subpath_contents = os.listdir(parent_directory + path_contents[cell_number])

            #  Figure out what cell number you are actually looking in...
            cell_folder_number = path_contents[cell_number]
            cell_folder_number = int(cell_folder_number[4:])

            file_number = 0
            for file_idx in subpath_contents:
                filename = subpath_contents[file_number]

                # Regular Expression to find out the channel.
                result_1 = re.search("1_CH0(\d)", filename)
                if result_1 is not None:
                    channel_number = result_1.group(1)

                    #  Regular Expression to find out the Time Point
                    result_2 = re.search("(\d\d\d\d\d\d)", filename)
                    if result_2 is not None:
                        time_number = result_2.group(1)

                        #  Regular Expression to Avoid Processing Data that Has Already Been Sheared
                        if not re.search("Shear", filename):

                            # Build up matrix for parallel processing
                            if loop_idx == 0:
                                input_arguments = np.array([cell_folder_number, channel_number, time_number])
                            else:
                                new_line = np.array([cell_folder_number, channel_number, time_number])
                                #  print(new_line)
                                input_arguments = np.concatenate((input_arguments, new_line), axis=0)
                            loop_idx += 1
                file_number += 1
        cell_number += 1
    size_of_arguments = int(np.size(input_arguments) / 3)
    parse_output = np.reshape(input_arguments, (size_of_arguments, 3))
    return parse_output

def deskew(inArray, angle, dz, xypixelsize):
    (z_len, y_len, x_len) = inArray.shape
    Trans = np.cos(angle * np.pi / 180) * dz / xypixelsize
    widenBy = np.uint16(np.ceil(z_len * np.cos(angle * np.pi / 180) * dz / xypixelsize))

    inArrayWiden = np.zeros((z_len, y_len, x_len + widenBy))
    inArrayWiden[:z_len, :y_len, :x_len] = inArray
    output = np.zeros((z_len, y_len, x_len + widenBy))

    xF, yF = np.meshgrid(np.arange(x_len + widenBy), np.arange(y_len))

    for k in range(z_len):
        inSlice = inArrayWiden[k, :, :]
        inSliceFFT = np.fft.fftshift(np.fft.fft2(inSlice))
        inSliceFFTTrans = inSliceFFT * np.exp(-1j * 2 * np.pi * xF * Trans * k / (x_len + widenBy))
        output_temp = np.abs(np.fft.ifft2(np.fft.ifftshift(inSliceFFTTrans)))
        output[k, :, :] = output_temp

    output[output < 0] = 0
    return np.uint16(output)  # return uint16 data to save as tiff

def process_image(image_info):
    cellidx, chidx, tidx = image_info
    chidx = int(chidx)
    tidx = int(tidx)

    imname = "1_CH0" + str(chidx) + "_" + str(("{:06d}".format(tidx))) + ".tif"
    imfile = parent_directory + "Cell" + str(cellidx) + "/" + imname

    export_name = imfile.replace('.tif', '') + '_fullShear.tif'
    exists = os.path.isfile(export_name)
    if exists:
        print("Shearing Complete Already")
    else:
        print(imfile)
        imarray = imread(imfile)
        imarray = deskew(imarray, angle, dz, xypixelsize)
        imsave(export_name, imarray)

input_arguments = parse_directory(parent_directory)
indices_to_delete = []
indexCounter = 0
for idx in input_arguments:
    temp = input_arguments[indexCounter]
    imname = "1_CH0" + str(int(temp[1])) + "_" + str(("{:06d}".format(int(temp[2])))) + ".tif"
    imfile = parent_directory + "Cell" + str(temp[0]) + "/" + imname
    export_name = imfile.replace('.tif', '') + '_fullShear.tif'
    exists = os.path.isfile(export_name)
    if exists:
        new_line = np.array([int(indexCounter)])
        indices_to_delete = np.concatenate((indices_to_delete, new_line), axis=0)
    indexCounter += 1

#  Delete the rows that already have been processed.
final_input_arguments = np.delete(input_arguments, indices_to_delete, 0)

if __name__ == '__main__':
    with Pool(numberThreads) as p:
        p.map(process_image, final_input_arguments)
print('Complete')
