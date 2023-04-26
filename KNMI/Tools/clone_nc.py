#!/usr/bin/env python
# http://opendap.knmi.nl/knmi/thredds/fileServer/radarprecipclim/RAD_NL21_RAC_MFBS_5min_NC/2003/04/RAD_NL21_RAC_MFBS_5min_200304302255.nc
import datetime as dt
import sys, os
import numpy as np
import netCDF4 as nc

inname = sys.argv[1]
basename,ext=os.path.splitext(inname)
outname = basename+'_max'+ext

xdimname = 'x'
ydimname = 'y'
timedimname = 'time'
outformat="NETCDF3_CLASSIC"      # explicit classic format 

IMIN = 0
IMAX = 1
IAVG = 2
ISTD = 3
IMED = 4
# outformat="NETCDF4"              # explicit NC4 format

#============================================================================

def find_unlimited_dimension(dsin):
    unlimdim = None
    unlimdimname = ''
    for dimname,dimobj in dsin.dimensions.iteritems():
        if (dimobj.isunlimited()):
            if (unlimdim is None):
               sys.stderr.write("Additional unlimited dimension [%s] !\n" % (dimname))
            else:
                unlimdim = dimobj
                unlimdimname = dimname
    return (dimname,dimobj)

def clone_minmax(dsin, dsout):
    # Copy dimensions
    global xdimname
    global ydimname
    global timedimname
    global outformat
    for dname, the_dim in dsin.dimensions.iteritems():
        dimlen = the_dim.size
        dsout.createDimension(dname, dimlen if not the_dim.isunlimited() else None)
    
    # Copy global attributes
    global_atts = ({k: dsin.getncattr(k) for k in dsin.ncattrs()})
    dsout.setncatts(global_atts)
    
    # Assume time dimension = unlimited dimension
    timedim, timedimobj = find_unlimited_dimension(dsin)

    # Copy variables
    for v_name, varin in dsin.variables.iteritems():
        v_name_out = v_name
        varout_datatype = varin.datatype
        if (outformat=='NETCDF3_CLASSIC'):
           if (varin.datatype=='uint16'): varout_datatype='int16'

        outVar = dsout.createVariable(v_name_out, varout_datatype, varin.dimensions)

        # Copy variable attributes
        ncatts = {}
        for attrib_name in varin.ncattrs():
           attrib_value = varin.getncattr(attrib_name)
           # Fillvalue and offset in the same type as the newly created datatype
           if (attrib_name in ['_FillValue','add_offset']):
              attrib_value = attrib_value.astype(outVar.datatype.name)
           ncatts[attrib_name] = attrib_value
        outVar.setncatts(ncatts)
        outVar[:] = varin[:]   # Only copy data for time-INdependent variables in this stage
    return
           
    
####################################################################################################

fnin = sys.argv[1]
fnout = sys.argv[2]
dsin=nc.Dataset(fnin,"r")
dsout=nc.Dataset(fnout,"w")
sys.stderr.write('Cloning file')
clone_minmax(dsin, dsout)

#   add variable stations
sys.stderr.write('Adding variable "stations"')
outVar = dsout.createVariable('stations', NF90_DOUBLE,(stations))
nstations=dat.dimensions['stations'].size
outVar[:] = np.array(range(nstations))+float(1)

dsout.close()
sys.stdout.write("Done!\n")

    
    



