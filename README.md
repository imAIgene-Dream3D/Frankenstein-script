# The issue

Imaris (Oxford instruments) is a widely used sofwtare for image analysis. recent version 10.10 released a pixel classifier to improve on object segmentation. However it only allows training the classifier on one file at a time, which limits a wide application of such classifier on different experiments.
# Frankenstein script

The Frankenstein script allows to train a pixel classifier of Imaris software (Oxford instruments) on more then one file.

This script takes timepoints from all the `.ims` files in a given directory, and generates a new file containing sections of every original file. 

There is no need to edit the script. All parameters can be given in the terminal or in an executable (using Pyinstaller).

## How to cite this script
If you find these extensions useful in your research, support our efforts by citing it as:
- E. Rios-Jimenez, "Frankenstein script", imAIgenelab 2024. https://www.imaigene-lab.com/

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

## Set-up & Requirements
The requirements are included in the `requirements.txt` file. You can run 
```{python}
pip install --no-cache-dir -r requirements.txt
```
on your environment to install the required packages.

This pipeline was run using python v 3.11.9. The installation of additional packages is required to run the script:
- `pyometiff` v 1.0.0
- `h5py` v 3.11.0
- `numpy` v 2.0.0


