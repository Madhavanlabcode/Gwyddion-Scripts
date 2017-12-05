### threeFoldSym.py
###
### This function allows quick three fold symmetrization of an image
### (intended to be of a FFT)
###
### USING PYGWY REQUIRES 32 BIT GWYDDION AND 32 BIT PYTHON 2.7
###
###

import gwy


plugin_menu="/Symmetrize/Three Fold Symmetrization"
plugin_type="PROCESS"

PI=3.14159265359

def run():
	### Create undo point
	#key=gwy.gwy_app_data_browser_get_current(gwy.APP_DATA_FIELD_KEY)
	#gwy.gwy_app_undo_checkpoint(gwy.data, [key])
	
	### get Data Field
	
	originalDField=gwy.gwy_app_data_browser_get_current(gwy.APP_DATA_FIELD)
	dummyDF=originalDField.new_alike()
	
	sixtyDegreeField=originalDField.new_rotated(dummyDF,60.*(PI/180.),gwy.INTERPOLATION_LINEAR,gwy.ROTATE_RESIZE_SAME_SIZE)
	
	onetwentyDegreeField=originalDField.new_rotated(dummyDF,120.*(PI/180.),gwy.INTERPOLATION_LINEAR,gwy.ROTATE_RESIZE_SAME_SIZE)
	
	originalDField.sum_fields(originalDField,sixtyDegreeField)
	originalDField.sum_fields(originalDField,onetwentyDegreeField)
	originalDField.multiply(1./3.)
	
	originalDField.data_changed()
