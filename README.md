# Frankenstein script

The Frankenstein script is intended to improve the pixel classifier training in Imaris software, from Oxford instruments


This script takes timepoints from all the `.ims` files in a given directory, and generates a new file containing sections of every original file. The intention is to better train the pixel classifier as to avoid overfitting.

There is no need to edit the script. All parameters can be given in the terminal or in an executable (using Pyinstaller).

## Input
The input will be the directory in which the `.ims` files are located

## Parameters
The user will need to introduce the following parameters:
- **Directory where the `.ims` files are located**. It will take all the .ims files in that directory)
- **Number of timepoints in the output file.** The user will be informed of the total number of timepoints in the processed `.ims` files and decide the final length
- **Voxel size (x, y, z)**. This information can be found in the properties of the original Imaris files (Ctrl + I/ Cmd + I in Imaris)

*Note that this script will only work if all the .ims files have the same number of channels and resolution (voxel size).*

## Output
The output is a file called `frankenstein_file.tif`, in 5D TIF format, that will be located in the provided input directory

To finish the conversion to .ims file, just upload it into the ImarisFileConverter and click `Start all`.

## Built-in functions
- `ims_to_numpy()`: Converts an `.ims` file into numpy array, based on the `ims_channels_to_h5()` function in general_segmentation_functions.image_handling
- `write_tiff()`: Converts a TIFF file from a numpy array. Provides names for the channels in Imaris (RS1, RS2, ...) and can receive voxel size dimension order information.
- `ims_files_in_wd()`: Provides a list with all the `.ims` files in a given directory
