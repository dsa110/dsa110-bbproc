import argparse
import numpy as np
import matplotlib.pyplot as plt
import sigproc
import warnings
from pathlib import Path
import sys

tstart = sys.argv[1];  # mjd trigger
beamID = sys.argv[2];  # beam index

filfile = sys.argv[3]; # spliced filterbank file

fnames=[sys.argv[3+1],\
sys.argv[3+2],\
sys.argv[3+3],\
sys.argv[3+4],\
sys.argv[3+5],\
sys.argv[3+6],\
sys.argv[3+7],\
sys.argv[3+8],\
sys.argv[3+9],\
sys.argv[3+10],\
sys.argv[3+11],\
sys.argv[3+12],\
sys.argv[3+13],\
sys.argv[3+14],\
sys.argv[3+15],\
sys.argv[3+16]];

nChans = 384;

fname = fnames[0];
print('processing file '+fname);
filsize = int(Path(fname).stat().st_size);
nSam = int(filsize/nChans);
data = np.fromfile(fname, dtype=np.uint8, count=-1, offset=0);
alldata = np.reshape(data,(nChans,nSam),order='F');

for fname in fnames[1:]:
    print('processing file '+fname);
    filsize = int(Path(fname).stat().st_size);
    nSam = int(filsize/nChans);
    data = np.fromfile(fname, dtype=np.uint8, count=-1, offset=0);
    data = np.reshape(data,(nChans,nSam),order='F');
    alldata = np.concatenate((alldata,data),axis=0);


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


    
fhead = {b'telescope_id': b'66',    # DSA?
  b'az_start': b'0.0',  # not sure where to find
  b'nbits': str(8).encode(),
  b'source_name': 'source'.encode(),
  b'data_type': b'1',   # ???
  b'nchans': str(384*16).encode(),
  b'machine_id': b'99', # ??
  b'tsamp': str(4*8.192e-6).encode(),
  b'foff': str(-0.03051757812).encode(),
  b'nbeams': str(beamID).encode(),
  b'fch1': str(1498.75).encode(),
  b'tstart': str(tstart).encode(),
  b'nifs': b'1'}

header_string = b''
header_string += to_sigproc_keyword(b'HEADER_START')
    
for keyword in fhead.keys():
    if keyword not in header_keyword_types.keys():
        pass
    else:
        header_string += to_sigproc_keyword(keyword, fhead[keyword])
    
header_string += to_sigproc_keyword(b'HEADER_END')

outfile = open(filfile, 'wb');
outfile.write(header_string);

np.reshape(alldata,(alldata.shape[0]*alldata.shape[1]),order='F').astype('uint8').tofile(outfile);

outfile.close();
