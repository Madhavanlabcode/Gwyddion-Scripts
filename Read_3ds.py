#### Read_3ds
### This is meant to read nanonis 3ds files to be used in
### the software Gwyddion. This is done through their pygwy module
### 
### USING PYGWY REQUIRES 32 BIT GWYDDION AND 32 BIT PYTHON 2.7
### 
### This program will utilize the nanonispy (made for Python 3.x) library by copying the necessary
### parts into this program (directory?) and changing them for Python 2.7



### Import necessary modules

import gwy #pygwy modules and functions

import os
import sys
import numpy as np

### CREATE INTERFACE TO LOAD load 3ds
#########################################################
_end_tags = dict(grid=':HEADER_END:', scan='SCANIT_END', spec='[DATA]')

class NanonisFile(object):

    """
    Base class for Nanonis data files (grid, scan, point spectroscopy).

    Handles methods and parsing tasks common to all Nanonis files.

    Parameters
    ----------
    fname : str
        Name of Nanonis file.

    Attributes
    ----------
    datadir : str
        Directory path for Nanonis file.
    basename : str
        Just the filename, no path.
    fname : str
        Full path of Nanonis file.
    filetype : str
        filetype corresponding to filename extension.
    byte_offset : int
        Size of header in bytes.
    header_raw : str
        Unproccessed header information.
    """

    def __init__(self, fname):
        self.datadir, self.basename = os.path.split(fname)
        self.fname = fname
        self.filetype = self._determine_filetype()
        self.byte_offset = self.start_byte()
        self.header_raw = self.read_raw_header(self.byte_offset)

    def _determine_filetype(self):
        """
        Check last three characters for appropriate file extension,
        raise error if not.

        Returns
        -------
        str
            Filetype name associated with extension.

        Raises
        ------
        UnhandledFileError
            If last three characters of filename are not one of '3ds',
            'sxm', or 'dat'.
        """

        if self.fname[-3:] == '3ds':
            return 'grid'
        elif self.fname[-3:] == 'sxm':
            return 'scan'
        elif self.fname[-3:] == 'dat':
            return 'spec'
        else:
            raise UnhandledFileError('{} is not a supported filetype or does not exist'.format(self.basename))

    def read_raw_header(self, byte_offset):
        """
        Return header as a raw string.

        Everything before the end tag is considered to be part of the header.
        the parsing will be done later by subclass methods.

        Parameters
        ----------
        byte_offset : int
            Size of header in bytes. Read up to this point in file.

        Returns
        -------
        str
            Contents of filename up to byte_offset as a decoded binary
            string.
        """

        with open(self.fname, 'rb') as f:
            return f.read(byte_offset).decode(encoding="latin1")

    def start_byte(self):
        """
        Find first byte after end tag signalling end of header info.

        Caveat, I believe this is the first byte after the end of the
        line that the end tag is found on, not strictly the first byte
        directly after the end tag is found. For example in Scan
        __init__, byte_offset is incremented by 4 to account for a
        'start' byte that is not actual data.

        Returns
        -------
        int
            Size of header in bytes.
        """

        with open(self.fname, 'rb') as f:
            tag = _end_tags[self.filetype]

            # Set to a default value to know if end_tag wasn't found
            byte_offset = -1

            for line in iter(f.readline,''):
                # Convert from bytes to str
                entry = line.strip().decode()
                if tag in entry:
                    byte_offset = f.tell()
                    break

            if byte_offset == -1:
                raise FileHeaderNotFoundError(
                        'Could not find the {} end tag in {}'.format(tag, self.basename)
                        )
        return byte_offset

class Grid(NanonisFile):

    """
    Nanonis grid file class.

    Contains data loading method specific to Nanonis grid file. Nanonis
    3ds files contain a header terminated by '\r\n:HEADER_END:\r\n'
    line, after which big endian encoded binary data starts. A grid is
    always recorded in an 'up' direction, and data is recorded
    sequentially starting from the first pixel. The number of bytes
    corresponding to a single pixel will depend on the experiment
    parameters. In general the size of one pixel will be a sum of

        - # fixed parameters
        - # experimental parameters
        - # sweep signal points (typically bias).

    Hence if there are 2 fixed parameters, 8 experimental parameters,
    and a 512 point bias sweep, a pixel will account 4 x (522) = 2088
    bytes of data. The class intuits this from header info and extracts
    the data for you and cuts it up into each channel, though normally
    this should be just the current.

    Currently cannot accept grids that are incomplete.

    Parameters
    ----------
    fname : str
        Filename for grid file.

    Attributes
    ----------
    header : dict
        Parsed 3ds header. Relevant fields are converted to float,
        otherwise most are string values.
    signals : dict
        Dict keys correspond to channel name, with values being the
        corresponding data array.

    Raises
    ------
    UnhandledFileError
        If fname does not have a '.3ds' extension.
    """

    def __init__(self, fname):
        _is_valid_file(fname, ext='3ds')
        super(Grid,self).__init__(fname)
        self.header = _parse_3ds_header(self.header_raw)
        self.signals = self._load_data()
        self.signals['sweep_signal'] = self._derive_sweep_signal()
        self.signals['topo'] = self._extract_topo()

    def _load_data(self):
        """
        Read binary data for Nanonis 3ds file.

        Returns
        -------
        dict
            Channel name keyed dict of 3d array.
        """
        # load grid params
        nx, ny = self.header['dim_px']
        num_sweep = self.header['num_sweep_signal']
        num_param = self.header['num_parameters']
        num_chan = self.header['num_channels']
        data_dict = dict()

        # open and seek to start of data
        f = open(self.fname, 'rb')
        f.seek(self.byte_offset)
        data_format = '>f4'
        griddata = np.fromfile(f, dtype=data_format)
        f.close()

        # pixel size in bytes
        exp_size_per_pix = num_param + num_sweep*num_chan
        
        # reshape from 1d to 3d
        griddata_shaped = griddata.reshape((nx, ny, exp_size_per_pix))

        # experimental parameters are first num_param of every pixel
        params = griddata_shaped[:, :, :num_param]
        data_dict['params'] = params

        # extract data for each channel
        for i, chann in enumerate(self.header['channels']):
            start_ind = num_param + i * num_sweep
            stop_ind = num_param + (i+1) * num_sweep
            data_dict[chann] = griddata_shaped[:, :, start_ind:stop_ind]

        return data_dict

    def _derive_sweep_signal(self):
        """
        Computer sweep signal.

        Based on start and stop points of sweep signal in header, and
        number of sweep signal points.

        Returns
        -------
        numpy.ndarray
            1d sweep signal, should be sample bias in most cases.
        """
        # find sweep signal start and end from a given pixel value
        sweep_start, sweep_end = self.signals['params'][0, 0, :2]
        num_sweep_signal = self.header['num_sweep_signal']

        return np.linspace(sweep_start, sweep_end, num_sweep_signal, dtype=np.float32)

    def _extract_topo(self):
        """
        Extract topographic map based on z-controller height at each
        pixel.

        The data is already extracted, though it lives in the signals
        dict under the key 'parameters'. Currently the 4th column is the
        Z (m) information at each pixel, should update this to be more
        general in case the fixed/experimental parameters are not the
        same for other Nanonis users.

        Returns
        -------
        numpy.ndarray
            Copy of already extracted data to be more easily accessible
            in signals dict.
        """
        return self.signals['params'][:, :, 4]

class UnhandledFileError(Exception):

    """
    To be raised when unknown file extension is passed.
    """
    pass


class FileHeaderNotFoundError(Exception):

    """
    To be raised when no header information could be determined.
    """
    pass


def _parse_3ds_header(header_raw):
    """
    Parse raw header string.

    Empirically done based on Nanonis header structure. See Grid
    docstring or Nanonis help documentation for more details.

    Parameters
    ----------
    header_raw : str
        Raw header string from read_raw_header() method.

    Returns
    -------
    dict
        Channel name keyed dict of 3d array.
    """
    # cleanup string and remove end tag as entry
    header_entries = header_raw.split('\r\n')
    header_entries = header_entries[:-2]

    # software version 'generic 5' had an extra header entry
    if header_entries[2] == 'Filetype=Linear':
        header_entries.pop(2)

    header_dict = dict()

    # grid dimensions in pixels
    dim_px_str = _split_header_entry(header_entries[0])
    header_dict['dim_px'] = [int(val) for val in dim_px_str.split(' x ')]

    # grid frame center position, size, angle
    grid_str = _split_header_entry(header_entries[1], multiple=True)
    header_dict['pos_xy'] = [float(val) for val in grid_str[:2]]
    header_dict['size_xy'] = [float(val) for val in grid_str[2:4]]
    header_dict['angle'] = float(grid_str[-1])

    # sweep signal
    header_dict['sweep_signal'] = _split_header_entry(header_entries[2])

    # fixed parameters
    header_dict['fixed_parameters'] = _split_header_entry(header_entries[3], multiple=True)

    # experimental parameters
    header_dict['experimental_parameters'] = _split_header_entry(header_entries[4], multiple=True)

    # number of parameters (each 4 bytes)
    header_dict['num_parameters'] = int(_split_header_entry(header_entries[5]))

    # experiment size in bytes
    header_dict['experiment_size'] = int(_split_header_entry(header_entries[6]))

    # number of points of sweep signal
    header_dict['num_sweep_signal'] = int(_split_header_entry(header_entries[7]))

    # channel names
    header_dict['channels'] = _split_header_entry(header_entries[8], multiple=True)
    header_dict['num_channels'] = len(header_dict['channels'])

    # measure delay
    header_dict['measure_delay'] = float(_split_header_entry(header_entries[9]))

    # metadata
    header_dict['experiment_name'] = _split_header_entry(header_entries[10])
    header_dict['start_time'] = _split_header_entry(header_entries[11])
    header_dict['end_time'] = _split_header_entry(header_entries[12])
    header_dict['user'] = _split_header_entry(header_entries[13])
    header_dict['comment'] = _split_header_entry(header_entries[14])

    return header_dict

def _split_header_entry(entry, multiple=False):
    """
    Split 3ds header entries by '=' character. If multiple values split
    those by ';' character.
    """

    _, val_str = entry.split("=")

    if multiple:
        return val_str.strip('"').split(';')
    else:
        return val_str.strip('"')
        
def save_array(file, arr, allow_pickle=True):
    """
    Wrapper to numpy.save method for arrays.

    The idea would be to use this to save a processed array for later
    use in a matplotlib figure generation scripts. See numpy.save
    documentation for details.

    Parameters
    ----------
    file : file or str
        File or filename to which the data is saved.  If file is a file-
        object, then the filename is unchanged.  If file is a string, a
        ``.npy`` extension will be appended to the file name if it does
        not already have one.
    arr : array_like
        Array data to be saved.
    allow_pickle : bool, optional
        Allow saving object arrays using Python pickles. Reasons for
        disallowing pickles include security (loading pickled data can
        execute arbitrary code) and portability (pickled objects may not
        be loadable on different Python installations, for example if
        the stored objects require libraries that are not available, and
        not all pickled data is compatible between Python 2 and Python
        3). Default: True
    """
    np.save(file, arr, allow_pickle=allow_pickle)
    
def load_array(file, allow_pickle=True):
    """
    Wrapper to numpy.load method for binary files.

    See numpy.load documentation for more details.

    Parameters
    ----------
    file : file or str
        The file to read. File-like objects must support the
    ``seek()`` and ``read()`` methods. Pickled files require that the
    file-like object support the ``readline()`` method as well.
    allow_pickle : bool, optional
        Allow loading pickled object arrays stored in npy files. Reasons
        for disallowing pickles include security, as loading pickled
        data can execute arbitrary code. If pickles are disallowed,
        loading object arrays will fail. Default: True

    Returns
    -------
    result : array, tuple, dict, etc.
        Data stored in the file. For ``.npz`` files, the returned
        instance of NpzFile class must be closed to avoid leaking file
        descriptors.
    """
    return np.load(file)
    
def _is_valid_file(fname, ext):
    """
    Detect if invalid file is being initialized by class.
    """
    if fname[-3:] != ext:
        raise UnhandledFileError('{} is not a {} file'.format(fname, ext))

### END INTERFACE TO LOAD 3ds
###############################################################



### NEED TO DEFINE THESE VARIABLES TO INTERFACE WITH GWYDDION

plugin_type = "FILE"
plugin_desc = "3ds file format used by Nanonis for Grid Spectroscopy"

### The pygwy interface requires four functions for loading and saving files
### They are detect_by_filename, detect_by_content, load, and save.
### Each function needs specific inputs and return types.

### Functions to determine if the filetype is valid and which load methods
### to use
def detect_by_filename(filename):
	if (filename.endswith(".3ds")):
		return 100
	else:
		return 0
		
def detect_by_content(filename, head, tail, filesize):
	if (head.startswith("Grid")):
		return 100
	else:
		return 0

### Load the file into the Gwyddion data types.		

def load(filename, mode=None):
	### Load returns a container object, initialize here
	mainC=gwy.Container()
	
	### Load the file into the Nanonispy Grid object
	grid=Grid(filename) 
	
	### Unpack the file for clarity, assuming the header format is constant.
	
	dim_x, dim_y=grid.header.get("dim_px")
	real_x, real_y= grid.header.get("size_xy")
	pos_x, pos_y = grid.header.get("pos_xy")
	channels=grid.header.get("channels")
	num_sweep_signal=grid.header.get("num_sweep_signal")
	num_parameters=grid.header.get("num_parameters")
	experimental_parameters=grid.header.get("experimental_parameters")
	sweep_size=abs(grid.signals.get("sweep_signal")[0]-grid.signals.get("sweep_signal")[-1])
	
	### Need to be able to load a arbitrary amount of 3D channel data
	for i in range(len(channels)):
		###Create the gwy Brick object and set dim/units
		tmpBrick=gwy.Brick(dim_x,dim_y,num_sweep_signal,real_x,real_y,sweep_size,1)
		tmpBrick.set_si_unit_x(gwy.SIUnit("m"))
		tmpBrick.set_si_unit_y(gwy.SIUnit("m"))
		tmpBrick.set_si_unit_z(gwy.SIUnit("V"))
		if "(V)" in channels[i]:
			tmpBrick.set_si_unit_w(gwy.SIUnit("V"))
		else:
			tmpBrick.set_si_unit_w(gwy.SIUnit("A"))
			
		### Transfer data to grid object	
		volData=grid.signals.get(channels[i])
		for x in range(dim_x):
			for y in range(dim_y):
				for z in range(num_sweep_signal):
					tmpBrick.set_val(x,y,z,volData[x,y,z])
					
		###Load Brick into main Container
		mainC.set_object_by_name("/brick/"+str(i),tmpBrick)
		
		###Create a preview by slicing the volume data
		tmpPreviewDataField=gwy.DataField(dim_x, dim_y, real_x, real_y,1)
		tmpBrick.extract_plane(tmpPreviewDataField,0,0,0,dim_x,dim_y,-1,True)
		mainC.set_object_by_name("/brick/"+str(i)+"/preview",tmpPreviewDataField)
		
		### set title of volume data
		mainC.set_string_by_name("/brick/"+str(i)+"/title",str(channels[i]))
		
		### Make sure data is visible upon loading 
		mainC.set_boolean_by_name("/brick/"+str(i)+"/visible",False)
		
		### TO DO: figure out how to store the meta data in Container
		
	###Load the topograph to display as well
	topoDataField=gwy.DataField(dim_x,dim_y,real_x,real_y,1)
	topo=grid.signals.get('topo')
	for x in range(dim_x):
		for y in range(dim_y):
			topoDataField.set_val(x,y,topo[x,y])
	
	topoDataField.set_si_unit_xy(gwy.SIUnit('m'))
	topoDataField.set_si_unit_z(gwy.SIUnit('m'))
	
	mainC.set_object_by_name("/0/data",topoDataField)
	mainC.set_string_by_name("/0/data/title","topo")
	mainC.set_boolean_by_name("/0/data/visible",True)
	
	return mainC
	
def save(data, filename, mode=None):
	return True
	
