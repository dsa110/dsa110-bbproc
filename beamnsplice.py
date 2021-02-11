# beamform and splices files together
# example:
# python beamnsplice.py 08feb21 /mnt/data/dsa110/T3/calibs/beamformer_weights_corrXX_06feb21.dat 264386032 4.813138888888889 9.45282 422.555 126 fl_2.out fl_2.out fl_2.out fl_2.out fl_2.out fl_2.out fl_2.out fl_2.out fl_2.out fl_2.out fl_2.out fl_2.out fl_2.out fl_2.out fl_2.out fl_2.out &

import os
import sys
import os.path
import glob
import numpy as np
import yaml
import astropy.units as u
from astropy.io import ascii
#import matplotlib.pyplot as plt
import json
import subprocess
from sigpyproc.Readers import FilReader
from pathlib import Path
import time

dirname  = sys.argv[1];
calfile  = sys.argv[2]; # contains XX instead of corr number
specnum  = int(sys.argv[3]);
timehr   = float(sys.argv[4]);
snr      = float(sys.argv[5]);
dm       = float(sys.argv[6]);
nBeamNum = int(sys.argv[7]);

# to test:
# inp = '261840032 4.7667777777777784 39.512 25.2644 136 fl_48.out fl_47.out fl_48.out fl_47.out fl_47.out fl_50.out fl_47.out fl_48.out fl_47.out fl_47.out fl_47.out fl_47.out fl_50.out fl_46.out fl_48.out fl_49.out';
# dirname = '09dec20';
# specnum = int(inp.split()[0]);
# dtime = float(inp.split()[1]);
# snr = float(inp.split()[2]);
# dm = float(inp.split()[3]);
# nBeamNum = int(inp.split()[4]);

nAnts = 24;     # antennas
nPols = 2;
nChans = 384;   # freq channels
freqs = [\
1498.75,\
1487.03125,\
1475.3125,\
1463.59375,\
1451.875,\
1440.15625,\
1428.4375,\
1416.71875,\
1405.0,\
1393.28125,\
1381.5625,\
1369.84375,\
1358.125,\
1346.40625,\
1334.6875,\
1322.96875];

# header stuff
header_keyword_types = {
    b'telescope_id' : b'<l',
    b'machine_id'   : b'<l',
    b'data_type'    : b'<l',
    b'barycentric'  : b'<l',
    b'pulsarcentric': b'<l',
    b'nbits'        : b'<l',
    b'nsamples'     : b'<l',
    b'nchans'       : b'<l',
    b'nifs'         : b'<l',
    b'nbeams'       : b'<l',
    b'ibeam'        : b'<l',
    b'rawdatafile'  : b'str',
    b'source_name'  : b'str',
    b'az_start'     : b'<d',
    b'za_start'     : b'<d',
    b'tstart'       : b'<d',
    b'tsamp'        : b'<d',
    b'fch1'         : b'<d',
    b'foff'         : b'<d',
    b'refdm'        : b'<d',
    b'period'       : b'<d',
    b'src_raj'      : b'<d',
    b'src_dej'      : b'<d',
#    b'src_raj'      : b'angle',
#    b'src_dej'      : b'angle',
    }

def to_sigproc_keyword(keyword, value=None):
    """ Generate a serialized string for a sigproc keyword:value pair
    If value=None, just the keyword will be written with no payload.
    Data type is inferred by keyword name (via a lookup table)
    Args:
        keyword (str): Keyword to write
        value (None, float, str, double or angle): value to write to file
    Returns:
        value_str (str): serialized string to write to file.
    """

    keyword = bytes(keyword)

    if value is None:
        return np.int32(len(keyword)).tostring() + keyword
    else:
        dtype = header_keyword_types[keyword]

        dtype_to_type = {b'<l'  : np.int32,
                         b'str' : str,
                         b'<d'  : np.float64}#,
#                         b'angle' : to_sigproc_angle}

        value_dtype = dtype_to_type[dtype]

        if value_dtype is str:
            return np.int32(len(keyword)).tostring() + keyword + np.int32(len(value)).tostring() + value
        else:
            return np.int32(len(keyword)).tostring() + keyword + value_dtype(value).tostring()


# processing
nFiles = 16;

fnames = [];
fcalib = [];
for k in range(nFiles):
    if k == 3:
        fnames.append('/home/user/data/T3/corr21/' + dirname + '/'+sys.argv[8+k]);
        fcalib.append(calfile.replace('XX',str(k+1).zfill(2)));
    else:
        fnames.append('/home/user/data/T3/corr' + str(k+1).zfill(2) + '/' + dirname + '/'+sys.argv[8+k]);
        fcalib.append(calfile.replace('XX',str(k+1).zfill(2)));

# measure size of the files, verify all files have the same size, estimate size of beamformed
FilSiz = np.zeros((nFiles));
for k in range(nFiles):
    FilSiz[k] = Path(fnames[k]).stat().st_size;
FilSiz /= (nAnts*nPols);

# run beamformer, wait until all files reached expected size

for k in range(nFiles):
    os.system('/home/user/T3_detect/code/beamformer -d ' + fnames[k] + ' -f ' + fcalib[k] + ' -o /home/user/data/findpulse/' + str(specnum) + '_corr' + str(k+1).zfill(2) + '_beam' + str(nBeamNum).zfill(3) + '.fil -z ' + str(freqs[k]) + ' -n ' + str(nBeamNum) + ' &');

time.sleep(5);

fullsize = 0.;
for k in range(nFiles):
    fullsize += Path('/home/user/data/findpulse/' + str(specnum) + '_corr' + str(k+1).zfill(2) + '_beam' + str(nBeamNum).zfill(3) + '.fil').stat().st_size;
    
while(fullsize < sum(FilSiz)):
    time.sleep(1);
    fullsize = 0.;
    for k in range(nFiles):
        fullsize += Path('/home/user/data/findpulse/' + str(specnum) + '_corr' + str(k+1).zfill(2) + '_beam' + str(nBeamNum).zfill(3) + '.fil').stat().st_size;

print('beamforming done, splicing...');

# splice files

fhead = {b'telescope_id': b'66',    # DSA?
  b'az_start': b'0.0',  # not sure where to find
  b'nbits': str(8).encode(),
  b'source_name': 'source'.encode(),
  b'data_type': b'1',   # ???
  b'nchans': str(nChans*nFiles).encode(),
  b'machine_id': b'99', # ??
  b'tsamp': str(4*8.192e-6).encode(),
  b'foff': str(-0.03051757812).encode(),
  b'nbeams': str(nBeamNum).encode(),
  b'fch1': str(freqs[0]).encode(),
  b'tstart': str(55000).encode(),
  b'refdm': str(dm).encode(),
  b'nifs': b'1'}

header_string = b''
header_string += to_sigproc_keyword(b'HEADER_START')
    
for keyword in fhead.keys():
    if keyword not in header_keyword_types.keys():
        pass
    else:
        header_string += to_sigproc_keyword(keyword, fhead[keyword])
    
header_string += to_sigproc_keyword(b'HEADER_END')

outfile = open('/home/user/data/findpulse/' + dirname + '_' + str(specnum) + '_beam' + str(nBeamNum).zfill(3) + '.fil', 'wb');
outfile.write(header_string);

alldata = [];
for k in range(nFiles):
    fname = '/home/user/data/findpulse/' + str(specnum) + '_corr' + str(k+1).zfill(2) + '_beam' + str(nBeamNum).zfill(3) + '.fil';
    data = np.fromfile(fname, dtype=np.uint8, count=-1, offset=0);
    data = np.reshape(data,(nChans,-1),order='F');
    alldata.append(data);


np.reshape(np.concatenate(np.asarray(alldata),axis=0),(-1),order='F').astype('uint8').tofile(outfile);

outfile.close();
os.system('rm /home/user/data/findpulse/*_corr*' + '_beam' + str(nBeamNum).zfill(3) + '.fil');

print('file spliced.')

# open file, plot spectrum  + time series

# nWin = 256;

# fname = '/home/user/T3_detect/processedcands/' + str(specnum) + '.fil';
# f = FilReader(fname);
# ts = f.dedisperse(dm);
# fts = ts.applyBoxcar(nWin);

# plt.figure()
# plt.plot(4*8.192e-6*np.arange(len(fts)),fts);
# plt.grid();
# plt.xlabel('time [s]');
# plt.ylabel('power [arbitrary]');
# plt.show();
