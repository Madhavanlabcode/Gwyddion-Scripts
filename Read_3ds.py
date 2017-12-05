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
import nanonispy27 as nap #used for reading 3ds files.

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
	grid=nap.read.Grid(filename) 
	
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
		mainC.set_boolean_by_name("/brick/"+str(i)+"/visible",True)
		
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
	
