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
The input will be the directory in which the `.ims` files are located. The python script will ask for the parameters during the execution.

In the Docker image, the parameters must be introduced when running the container. 

## Parameters
The user will need to introduce the following parameters:
- `working_directory` **Directory where the `.ims` files are located**. It will take all the .ims files in that directory)
- `size` **Number of timepoints in the output file.** The user will be informed of the total number of timepoints in the processed `.ims` files and decide the final length
- `x_vox, y_vox, z_vox` **Voxel size (x, y, z)**. This information can be found in the properties of the original Imaris files (Ctrl + I/ Cmd + I in Imaris)

*Note that this script will only work if all the .ims files have the same number of channels and resolution (voxel size).*

## Output
The output is a file called `frankenstein_file.tif`, in 5D TIF format, that will be located in the provided input directory

To finish the conversion to .ims file, just upload it into the ImarisFileConverter and click `Start all`.

## Docker image
A docker image of the script is provided, with all the requirements met and ready to use. You can find it as **erios12/frankenstein_script** [Click here to acces Docker hub](https://hub.docker.com/layers/erios12/frankenstein_script/1.0.0/images/sha256:2e3a3118d5422f7345202dbc1d4956b9b0e8345fb2325e6d3fdf8d4f6c06591d?uuid=F0C28C93-5A7A-4EB8-AFF3-B6DEB8B568A6)

The use of this docker image requires the parameters to be set **before** the execution, in the `docker run` command. You can copy, adapt and paste this command into your CMD/Powershell/WSL/console once the docker image has been built. You must have the docker environment and paths installed.

**Example:**
```
docker run --rm -v /c/Users/User/frankenstein_script/data:/app/data erios12/frankenstein_script:1.0.0 --size 200 --x_vox 0.5 --y_vox 0.5 --z_vox 1.0
```
- Parameters:
  - `working_directory`: Provided through the `-v` option. Basically copies the working directory into the `/app/data` directory in the docker Virtual Machine, so you can work with tour data in Docker images.
  - `size`: provided using the `--size` parameter option
  - `voxel_size`: provided using the `--x_vox`, `y_vox`, and `z_vox` parameter options (0.0 format)
- `-rm` option is to delete the container once the execution is completed

## Walkthrough video
We provide a walkthrough video for both the python script and the Docker implementation, to help resolve any issues. If you are having trouble, check the [Troubleshooting](#Troubleshooting) section.

https://github.com/user-attachments/assets/8c05c645-8844-4579-88bd-821f1cf53289


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

## Troubleshooting
- Make sure that all the `.ims` files in the directory are to be used to create the frankenstein file. They all must have the same resolution and channels.
- Do not insert a `size` value that exceeds the total sum of timepoints in the files provided.
- Make sure your machine has enough RAM and CPU power to work with image files.
- In the docker image implementation, in Windows, you can modify the RAM and CPU usage by creating a `.wslconfig` file. [More info](https://learn.microsoft.com/en-us/windows/wsl/wsl-config#wslconfig)
- Due to be working in a virtual environment, the Docker image can have a worse performance than running the raw script in your machine
- If the Docker container stops working without any errors, it has most likely run out of RAM memory.

