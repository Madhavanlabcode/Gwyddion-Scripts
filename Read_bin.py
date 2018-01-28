### Read_LayerBin
###
### This is meant to read layer bin files from the java code for Gwyddion
### 
### USING PYGWY REQUIRES 32 BIT GWYDDION AND 32 BIT PYTHON 2.7
###
### 


### Import necessary modules

import gwy #pygwy modules and functions
import numpy as np #easy storing of data
import struct #reading binary file


plugin_type = "FILE"
plugin_desc = "Layer bin file from Java Code"


def detect_by_filename(filename):
	if (filename.endswith('.bin')):
		return 100
	else:
		return 0
		
def detect_by_content(filname, head, tail, filesize):
	nx=int(struct.unpack('>i',head[0:4])[0])
	ny=int(struct.unpack('>i',head[4:8])[0])
	if (filesize==(24+nx*ny*8+nx*8+ny*8)):
		return 100
	else:
		return 0


def load(filename):
	### Load returns a container object, initialize here
	mainC=gwy.Container()
	
	
	with open(filename,"rb") as f:
		
		###Read the 'header' parts of the bin file
		nx=int(struct.unpack('>i',f.read(4))[0])
		ny=int(struct.unpack('>i',f.read(4))[0])
		bias=struct.unpack('>d',f.read(8))[0]
		setCurrent=struct.unpack('>d',f.read(8))[0]
		
		#load the metadata
		
		metaC=gwy.Container()
		
		
		metaC.set_string_by_name("Bias: ", str(bias))
		mainC.set_object_by_name("/0/meta",metaC)
		
		#initialize the dimensions and data arrays
		data=np.zeros((nx,ny))
		xspacing=np.zeros(nx)
		yspacing=np.zeros(ny)
		
		for x in range(nx):
			xspacing[x]=struct.unpack('>d',f.read(8))[0]
			
		for y in range(ny):
			yspacing[y]=struct.unpack('>d',f.read(8))[0]
		
		for x in range(nx):
			for y in range(ny):
				data[x,y]=struct.unpack('>d',f.read(8))[0]
				
		dField=gwy.DataField(nx,ny,abs(xspacing[0]-xspacing[-1]),abs(yspacing[0]-yspacing[-1]),True)
		dField.set_si_unit_xy(gwy.SIUnit('m'))
		dField.set_si_unit_z(gwy.SIUnit('V'))
		
		for x in range(nx):
			for y in range(ny):
				dField.set_val(x,y,data[x,y])
				
		mainC.set_object_by_name("/0/data",dField)
		mainC.set_boolean_by_name("/0/data/visible",True)
		
		return mainC
		
def save(data, filename, mode=None):
	return True		
				
