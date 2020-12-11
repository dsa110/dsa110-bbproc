# dsa110-bbproc

contains scripts and code for baseband processing

# beamformer.c
produces beamformed binary file from voltage data
	   -d voltage data file name [no default]
	   -f calibration file name[no default]
	   -o output file name[no default]
	   -z fch1 in MHz [default 1530]
	   -s interbeam separation in arcmin [default 1.4]
	   -n beam number [0 -- 255, default 127]
	   -h print usage

# findfile.py
finds a voltage file associated with a trigger (T2 specnum)
 [1] directory with all voltage files and dump.dat
 [2] T2 specnum associated with trigger of interest
