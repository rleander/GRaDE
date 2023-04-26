#!/usr/bin/env python
import os
import sys
import argparse as ap
import pprint as pp
import subprocess as sp         # redundant 

# TODO:
# Generic script that extracts a block of blocksize years of a combined (multi-site) .pr, .tg and .ev file  
# into temporary files and performs an action (in this case creating a netCDF) 
# by launching a separate program (in this case ncGRaDE), as a subprocess

#----------------------------------------------------------------------------------------------------
# Reads a file by blocks of 1000 year
# Runs the transformation program for rr and tm
# puts the blocks back together
programname='ncGRaDE0'
basename='simdata'
basename_out='simdata'
offsetyr=2000               # shift year by ...
refyr=2001                  # reference year of time in udunitstring (years since ...)
catchments='catchments.txt' # list of catchments
metafile='meta.txt'         # meta data file, global attributes
blocksize = 2000

# Command-line 
#   ./ncGRaDE -sim ${simbasename}_0${iblock} -start ${firstyr}0101 \
#                                            -end ${lastyr}1231    \
#                                            -offset ${offsetyr}   \
#                                            -refyear ${refyr}     \
#                                            -list ${catchments}   \
#                                            -meta ${metafile} -nc ${simbasename}_0${iblock}
#
#----------------------------------------------------------------------------------------------------


rrfilename = basename+'.rr'
frr=open(rrfilename,"r")
tgfilename = basename+'.tg'
ftg=open(tgfilename,"r")
evfilename = basename+'.ev'
fev=open(evfilename,"r")

blk = -1
lineno_rr=0
lineno_tg=0
lineno_ev=0
try:
    while True: 
        while True:
           line_rr = frr.readline().strip()
           lineno_rr = lineno_rr + 1
           if (line_rr[0]!='#'): break 
        while True:
           line_tg = ftg.readline().strip()
           lineno_tg = lineno_tg + 1
           if (line_rr[0]!='#'): break 
        while True:
           line_ev = fev.readline().strip()
           lineno_ev = lineno_ev + 1
           if (line_rr[0]!='#'): break 
        blk0 = blk
        date    = int(line_rr.split()[0])
        date_tg = int(line_tg.split()[0])
        date_ev = int(line_ev.split()[0])
        if ((date_tg!=date) or (date_ev!=date)):  # date inconsistency 
            sys.stderr.write("Date inconsistency!\n")
            sys.stderr.write("  line %5d : %16d in %s\n" % (lineno_rr,date,rrfilename))
            sys.stderr.write("  line %5d : %16d in %s\n" % (lineno_tg,date_tg,tgfilename))
            sys.stderr.write("  line %5d : %16d in %s\n" % (lineno_ev,date_ev,evfilename))
            sys.exit()
        year = int(date/10000)
        blk = int((year-1)/blocksize)
        if (blk!=blk0):
            if (blk0>=0):                         # valid block finished, then call transformation program
                ftmp_rr.close()                   # stop writing temporary file
                ftmp_tg.close()                   # stop writing temporary file
                ftmp_ev.close()                   # stop writing temporary file
                lastyear=year-1
                lastdate="%d1231"%lastyear
                # Build command-line
                transcmdl = programname+" "+" -sim tmpseries " \
                          + "-start %s "%firstdate  \
                          + "-end %s "%lastdate     \
                          + "-offset %d "%offsetyr  \
                          + "-refyear %d "%refyr    \
                          + "-list %s "%catchments  \
                          + "-meta %s "%metafile    \
                          + "-nc %s_%03d.nc "%((basename_out),blk0) \
                          + " ".join(sys.argv[4:])
                print transcmdl
                # Execute command-line
    #           sys.stderr.write("processing block %d: %s\n" % (blk0,transcmdl))        ! RL666
    #           os.system(transcmdl)                                                    ! RL666
            ftmp_rr = open('tmpseries.pr',"w")     # (re-)open the temporary file, overwrite mode
            ftmp_tg = open('tmpseries.tg',"w")     # (re-)open the temporary file, overwrite mode
            ftmp_ev = open('tmpseries.ev',"w")     # (re-)open the temporary file, overwrite mode
            firstyear=year
            firstdate="%d0101"%firstyear
        ftmp_rr.write("%s\n"%line_rr)              # copy lines
        ftmp_tg.write("%s\n"%line_tg)              # copy lines
        ftmp_ev.write("%s\n"%line_ev)              # copy lines
except:
    pass
ftmp_rr.close()
ftmp_tg.close()
ftmp_ev.close()
# Call transformation program again, last time
transcmdl = programname+" "+" -sim tmpseries " \
          + "-start %s "%firstdate  \
          + "-end %s "%lastdate     \
          + "-offset %d "%offsetyr  \
          + "-refyear %d "%refyr    \
          + "-list %s "%catchments  \
          + "-meta %s "%metafile    \
          + "-nc %s_%03d.nc "%((basename_out),blk0) \
          + " ".join(sys.argv[4:])
print transcmdl
os.system(transcmdl)
        
fin.close()
