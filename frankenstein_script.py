import sys
import os
import glob
from random import randint
from time import sleep
import numpy as np
import pyometiff as ometiff

class Image(object):
    def __init__(self, path, permission='r'):
        import os
        import sys
        path=str(path)
        self.known_image_types=[".ims"]
        self.path = os.path.realpath(path).replace("\\", "/")
        self.format = self.get_format() or ""
        self.permission = permission
        self.file=None
        self.ds=None
        self.image=None       
        self.dims=None
        self.axlab=["z", "y", "x"]
        self.elsize=[1,1,1]
        self.rel_elsize=[1,1,1]
        self.metadata={}
        self.pathparts=self.split_path()
        self.h5path=self.pathparts
        
        self.ds=None
        self.image=None
        
        if not os.path.exists(self.pathparts['file']) and self.permission=="r":
            print("WARNING: Path to image file does not exist\n{}".format(self.pathparts['file']))
            sys.exit()
    def ims_remove_zeropadding(self, in_place=False):
        unzeropadded_X=self.attr2datatype(self.file['/DataSetInfo/Image'].attrs['X'], int)
        unzeropadded_Y=self.attr2datatype(self.file['/DataSetInfo/Image'].attrs['Y'], int)
        unzeropadded_Z=self.attr2datatype(self.file['/DataSetInfo/Image'].attrs['Z'], int)
        try:
            self.image[:]
        except:
            print("No dataset has been loaded yet")
        else:
            if in_place:
                resize_list = tuple(None) * (self.image[:].ndim-3) + (unzeropadded_Z, unzeropadded_Y,unzeropadded_X)
                self.ds.resize(resize_list)
                self.image=self.ds[:]
                self.dims=self.ds.shape
            else:
                resize_list=(slice(None),) * (self.image[:].ndim-3) + (slice(0,unzeropadded_Z), slice(0,unzeropadded_Y), slice(0,unzeropadded_X))
                self.image=self.image[resize_list]

    def close(self):
        """Close a file."""
        import h5py
        try:
            if isinstance(self.file, h5py.File):
                self.file.close()
        except:
            print("Can't close file")

    def attrs2dict(self, attrs):
        dict={}
        for k,v in attrs.items():
            dict[k]=self.attr2str(v)
        return(dict)
    
    def attr2datatype(self, attrs,datatype=str):
        fixed_attrs=[]
        try:
            for t in attrs:
                try:
                    t=t.decode('utf-8')
                except:
                    pass
                finally:
                    fixed_attrs.append(t)

            fixed_attrs=datatype(''.join(fixed_attrs))
        except:
            fixed_attrs=attrs
        return(fixed_attrs)
    
    def attr2str(self, attrs):
        fixed_attrs=[]
        try:
            for t in attrs:
                try:
                    t=t.decode('utf-8')
                except:
                    pass
                finally:
                    fixed_attrs.append(t)
            fixed_attrs=''.join(fixed_attrs)
        except:
            fixed_attrs=attrs
        return(fixed_attrs)
    
    def get_axlab(self):
        """Get the axis labels."""

        if ((self.axlab is not None) and (len(self.axlab) > 0)):
            axlab = ''.join(self.axlab)
#             print("""WARNING:
#                   axlab already specified as {}""".format(self.axlab))
            return axlab

        formats = {
            '.h5': self.h5_get_axlab,
            '.ims': self.ims_get_axlab,
        }
        axlab = formats[self.format]()

        if axlab is None:
            axlab = 'zyxct'[:self.get_ndim()]
#             raise Exception("""WARNING: axlab is None;
#                                replaced by {}""".format(axlab))

        return axlab
    
    def get_format(self):
        import re

        image_type=None
        for ext in self.known_image_types:
            if re.match(".+"+ext+"/*.*$", self.path):
                image_type=ext
                break

        if not image_type:
            print(".{} is not a known image type".format(".".join(self.path.split(".")[1:])))
            print("Known image types are:", self.known_image_types)
        else:
            return(image_type)

    def split_path(self):
        file_comps={}
        file_comps['ext']=self.format
        try:
            before_ext, after_ext = self.path.split(self.format)
        except:
            before_ext = self.path
            after_ext=""
        file_comps['file'] = before_ext + file_comps['ext']
        file_comps['file_no_ext'] = before_ext
        if file_comps['ext']=="":
            file_comps['folder']=before_ext
            file_comps['filename']=""
            file_comps['filename_no_ext']=""
        else:
            file_comps['folder'], file_comps['filename'] = file_comps['file'].rsplit("/", 1)
            file_comps['filename_no_ext']=file_comps['filename'].split(file_comps['ext'])[0]
        file_comps['folder']=file_comps['folder']+"/"
        file_comps['substruct']=None
        if after_ext!='' and after_ext!='/':
            file_comps['substruct']=after_ext
            file_comps['group'], file_comps['dataset'] = after_ext.rsplit("/", 1)
        return(file_comps)
    
    def load(self, substruct='', channel=''):
        import os
        self.ds=None
        self.image=None

        image_loaders={
        '.ims' : self.ims_load
        }

        try:
            loader=image_loaders[self.format]
        except:
            print("No loader for {} exists yet".format(self.format))
        else:
            loader(substruct=substruct, channel=channel)

        self.axlab=self.get_axlab()
        return(self.image)

    def ims_get_elsize(self):
        im_info = self.file['/DataSetInfo/Image']

        extmin0 = float(self.attr2str(im_info.attrs['ExtMin0']))
        extmin1 = float(self.attr2str(im_info.attrs['ExtMin1']))
        extmin2 = float(self.attr2str(im_info.attrs['ExtMin2']))
        extmax0 = float(self.attr2str(im_info.attrs['ExtMax0']))
        extmax1 = float(self.attr2str(im_info.attrs['ExtMax1']))
        extmax2 = float(self.attr2str(im_info.attrs['ExtMax2']))

        extX = extmax0 - extmin0
        extY = extmax1 - extmin1
        extZ = extmax2 - extmin2

        dims = self.dims or self.ims_get_dims()

        elsizeX = extX / dims[2]
        elsizeY = extY / dims[1]
        elsizeZ = extZ / dims[0]

        elsize=[elsizeZ, elsizeY, elsizeX]
        self.elsize=elsize
        return(elsize)
    
    def get_rel_elsize(self, elsize):
        min_size=min(elsize)
        elsize_scaled= [round(x/min_size, 2) for x in elsize]
        return(elsize_scaled)
    
    def ims_load(self, remove_zeropadding=True, channel = '', substruct = ''):
        import h5py
        import numpy as np
        channel=str(channel)
        if self.file is None:
            self.file = h5py.File(name=self.pathparts['file'], mode=self.permission)

        substruct = substruct or self.pathparts['substruct']

        if substruct:
            self.pathparts['substruct']=substruct
            self.pathparts['group'], self.pathparts['dataset'] = self.pathparts['substruct'].rsplit("/", 1)
            try:
                self.ds = self.file[substruct]
            except:
                print("{} is not a valid channel substructure for the .ims image".format(substruct))
            else:
                self.metadata.update(self.attrs2dict(self.file[self.pathparts['group']].attrs))
                self.image=self.ds[:]
                self.dims = self.ims_get_dims()
                self.elsize=self.ims_get_elsize()
                self.rel_elsize=self.get_rel_elsize(self.elsize)
                self.dtype = self.ds.dtype
                self.chunks = self.ds.chunks
                if self.image.shape != self.dims[-(self.image.ndim):]:
                    print("Imaris image is  zeropadded: org({}), unzeropadded({})".format(self.image.shape, self.dims[:3]))
        elif channel:
            default_channels_path='/DataSet/ResolutionLevel 0/'
            
            try:
                nr_timepoints = len([x for x in self.file[default_channels_path]])
                if nr_timepoints>1:
                    images = []
                    print(f"Loading {nr_timepoints} timepoints...")
                    for frame in range(nr_timepoints):
                        self.pathparts['substruct']="{}/TimePoint {}/Channel {}/Data".format(default_channels_path, frame, channel)
                        self.pathparts['group'], self.pathparts['dataset'] = self.pathparts['substruct'].rsplit("/", 1)
                        # print("Loading TimePoint {}".format(frame))
                        self.ds=self.file[ self.pathparts['substruct']]
                        images.append(self.file[ self.pathparts['substruct']][:])
                    self.image = np.stack(images, axis=0)
                else:
                    channel_path = default_channels_path+"TimePoint 0/"
                    self.pathparts['substruct']="{}/Channel {}/Data".format(channel_path, channel)
                    self.pathparts['group'], self.pathparts['dataset'] = self.pathparts['substruct'].rsplit("/", 1)
                    self.ds = self.file[self.pathparts['substruct']]
                    self.image=self.ds[:]
            except Exception as error:
                print('An exception occurred: {}'.format(error))
                print("{} is not a valid channel within the .ims image".format(channel))
            else:
                self.metadata.update(self.attrs2dict(self.file[self.pathparts['group']].attrs))
                
                self.dims = self.ims_get_dims()
                self.elsize=self.ims_get_elsize()
                self.rel_elsize=self.get_rel_elsize(self.elsize)
                self.dtype = self.ds.dtype
                self.chunks = self.ds.chunks
                if self.image.shape != self.dims[-(self.image.ndim):]:
                    print("Imaris image is zeropadded: org({}), unzeropadded({})".format(self.image.shape, self.dims[-(self.image.ndim):]))
        else:
            print("Base .ims file loaded as no substructure is specified")

        if remove_zeropadding:
            self.ims_remove_zeropadding()
        return(self.image)

    
    def ims_list_channels(self):
        default_channels_path='/DataSet/ResolutionLevel 0/TimePoint 0/'
        default_info_path='/DataSetInfo'
        print(self.file[default_channels_path])
        channel_names={}
        for channel in self.file[default_channels_path]:
            try:
                channel_name=self.attr2str(self.file[default_info_path][channel].attrs["Name"])
                channel_names[channel.split(" ")[1]]=channel_name
                print("{} - {}".format(channel, channel_name))
            except:
                print(channel)
        return(channel_names)

    def ims_get_dims(self):

        if (self.ds is not None and self.metadata):
            #dimX=int("".join([str(x.decode()) for x in list(self.metadata["ImageSizeX"])]))
            dimX=self.attr2datatype(self.metadata["ImageSizeX"], int)
            dimY=self.attr2datatype(self.metadata["ImageSizeY"], int)
            dimZ=self.attr2datatype(self.metadata["ImageSizeZ"], int)

        else:
            dimX=self.attr2datatype(self.file['/DataSetInfo/Image'].attrs['X'], int)
            dimY=self.attr2datatype(self.file['/DataSetInfo/Image'].attrs['Y'], int)
            dimZ=self.attr2datatype(self.file['/DataSetInfo/Image'].attrs['Z'], int)
        dimC = len(self.file['/DataSet/ResolutionLevel 0/TimePoint 0/'])
        dimT = len(self.file['/DataSet/ResolutionLevel 0/'])

        dims=[dimT, dimZ, dimY, dimX]
        return(dims)
        

##### END OF Image CLASS ######

def get_image(filepath, substruct="", channel=""):
    channel=str(channel)
    if substruct:
        im = Image('{}/{}'.format(filepath, substruct), permission='r')
    else:
        im = Image(filepath, permission="r")
    im.load(channel=channel)
    if im.image is None:
        sys.exit("Image was not loaded correctly. Exiting...")
    if im.format==".ims":
        im.ims_remove_zeropadding(in_place=False)
    data = im.image
    im.close()
    return data

def ims_list_channels(self):
        default_channels_path='/DataSet/ResolutionLevel 0/TimePoint 0/'
        default_info_path='/DataSetInfo'
        print(self.file[default_channels_path])
        channel_names={}
        for channel in self.file[default_channels_path]:
            try:
                channel_name=self.attr2str(self.file[default_info_path][channel].attrs["Name"])
                channel_names[channel.split(" ")[1]]=channel_name
                print("{} - {}".format(channel, channel_name))
            except:
                print(channel)
        return(channel_names)

def ims_to_numpy(
        ims_path,
        h5_group_names=[]
        ):

    from frankenstein_script import Image, get_image
    import os
    from pathlib import Path
    ims_path=Path(ims_path)
    if os.path.exists(ims_path):
        image=Image(ims_path, 'r')
        image.load()
        ims_channels=image.ims_list_channels()
        if not h5_group_names:
            h5_group_names=["raw_channels"] * len(ims_channels)
        elif isinstance(h5_group_names, str):
            h5_group_names=[h5_group_names] * len(ims_channels)
        elif isinstance(h5_group_names, list):
            pass
        else:
            print("h5_group_names of unknown type: {}".format(type(h5_group_names)))

        print(ims_channels)

        # Extract channel information
        extract_zip=zip(ims_channels.keys(), ims_channels.values(), h5_group_names)
        raw_ch = []
        for channel in extract_zip:
            idx=channel[0]
            name=channel[1]
            group_name=channel[2]
            print("Saving channel {}".format(str(idx)))
            ims_image = get_image(ims_path, channel=str(idx))
            group=group_name
            substruct="{}/{}".format(group, name)
            raw_ch.append(ims_image)
        raw_ch = np.array(raw_ch)
        print(raw_ch.shape)
        image.close()

    else:
        print(f"Provided path does not exist:\n{ims_path}")

    return len(raw_ch[0]), raw_ch


def write_tiff(out_path, frank_data, voxel_size=['', '', ''], dim_order = "TZCYX"):
    
    # Metadata information
    metadata = {
        "PhysicalSizeX" : f"{voxel_size[0]}",
        "PhysicalSizeXUnit" : "µm",
        "PhysicalSizeY" : f"{voxel_size[1]}",
        "PhysicalSizeYUnit" : "µm",
        "PhysicalSizeZ" : f"{voxel_size[2]}",
        "PhysicalSizeZUnit" : "µm",
        "Channels" : {}
    }

    # Channel Color and name
    for ch in range(1, len(frank_data[0][0])+1):
        metadata['Channels'][f'{ch}'] = {}
        metadata['Channels'][f'{ch}']["Name"] =  f'RS{ch}'
        if ch == 1:
            metadata['Channels'][f'{ch}']["Color"] =  "(0, 255, 255)" # Cyan
        elif ch == 2:
            metadata['Channels'][f'{ch}']["Color"] =  "(255, 255, 0)" # Yellow
        elif ch == 3:
            metadata['Channels'][f'{ch}']["Color"] =  "(255, 0, 0)" # Red
        else:
            metadata['Channels'][f'{ch}']["Color"] =  "(255, 255, 255)" # White

    # Create an OME-TIFF writer
    writer = ometiff.OMETIFFWriter(
        fpath=out_path,
        dimension_order=dim_order,
        array=frank_data,
        metadata=metadata,
        explicit_tiffdata=False,
        bigtiff=True)

    writer.write()

def ims_files_in_wd(working_directory):

    # Verify if the entered directory is valid
    if not os.path.isdir(working_directory):
        print(f"The directory '{working_directory}' does not exist.")
        sys.exit(1)

    # Find all .ims files in the current directory
    ims_files = glob.glob(os.path.join(working_directory, '*.ims'))
    print(f"Found .ims files: {ims_files}")

    if len(ims_files) == 0:
        print("No '.ims' files found. Please choose antother directory")
        sleep(5)
        sys.exit()

    return ims_files

def main():
    # Get the working directory from the user
    working_directory = input("Please enter the working directory (copy from window explorer bar): ")
    working_directory = working_directory.replace("\\", "/")

    # Get ims_files from wd
    ims_files = ims_files_in_wd(working_directory)

    # Import ims files to numpy array/list
    data = []
    time = 0
    for ims_file in ims_files:
        tp, raw_ch = ims_to_numpy(ims_file)
        time += tp
        ch_reord = np.moveaxis(raw_ch, [0, 1], [1, 0])
        data.append(ch_reord)
    
    print(f"Total number of timepoints: {time}")

    # Frankenstein file parameters
    size = input("Please select the number of timepoints desired in your frankenstein file (200 recommended):")
    n_timepoints = int(size)//len(ims_files)
    if n_timepoints < 10:
        n_timepoints = 10
    print(f"{n_timepoints} timepoints will be randomly selected from every file")
    print(f"The resulting file will have {n_timepoints*len(ims_files)} timepoints ({n_timepoints} from each of the {len(ims_files)} '.ims' files)")

    # Generate final dataset
    filter_data = []
    for file in data:
        start = randint(0, len(file)- (n_timepoints+1))
        end = start + n_timepoints
        filter_data.append(file[start:end])
    
    # Concatenate along the first axis (timepoints)
    concat_data = np.concatenate(filter_data, axis=0)
    print("(T, C, Z, Y, X)")
    print(concat_data.shape)

    # Reorder shape for input in ImarisFileConverter
    print("Original shape (T, C, Z, Y, X)")
    frank_data = np.moveaxis(concat_data, [0, 1, 2], [0, 2, 1])
    print(frank_data.shape)
    print("Reordered shape for Imaris import (T, Z, C, Y, X)")

    
    # Voxel size
    x_vox = input("X voxel size (use 0.0 format):")
    y_vox = input("Y voxel size (use 0.0 format):")
    z_vox = input("Z voxel size (use 0.0 format):")
    voxel_size = [x_vox, y_vox, z_vox]

    print(f"Saving new TIFF channels in {working_directory}")

    # Save the NumPy array as a TIFF file
    out_path = os.path.join(working_directory, "frankenstein_file.tif")

    # Incorporates voxel_size and auto names channels
    write_tiff(out_path, frank_data, voxel_size)

    print("New frankenstein file in TIFF format completed. Use ImarisFileConverter to finish the connversion to .ims")
    print("Open the TIF file in the ImarisFileCOnverter and click 'Start all'")
    sleep(10)

# Execute only if main program (not as module)
if __name__ == "__main__":
    main()