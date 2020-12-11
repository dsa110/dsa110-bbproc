import numpy as np
import numpy.matlib
import sigproc
import matplotlib.pyplot as plt
from pathlib import Path
import sys

fname = sys.argv[1];    #file name

nChans = 384*16;

header_len = sigproc.len_header(fname);
filsize = int(Path(fname).stat().st_size);
nSam = int((filsize - header_len)/nChans);
data = np.fromfile(fname, dtype=np.uint8, count=-1, offset=header_len);

alldata = np.reshape(data,(nChans,nSam),order='F');
#alldata = alldata.astype(np.float64);

# medspec = np.median(alldata,axis=1);

# for k in range(nSam):
    # alldata[:,k] /= medspec;

# varspec = np.var(alldata,axis=1);

# for k in range(nSam):
    # alldata[:,k] /= varspec;

dataint = np.zeros((32,2048));
for k in range(32):
	tmp = np.mean(alldata[k*192:(k+1)*192,:],axis=0);
	for kk in range(2048):
		dataint[k,kk] = np.mean(tmp[kk*60:(kk+1)*60]);

plt.figure();
plt.subplot(311);
plt.imshow(dataint,aspect='auto',extent=[0,4,1487.03125,1322.9875-(1334.6875-1322.96875)]);
plt.xlabel('time [s]');
plt.ylabel('frequency [MHz]');
plt.subplot(312);
plt.plot(np.linspace(0,4,2048),np.mean(dataint,axis=0));
plt.xlabel('time [s]');
plt.ylabel('power');
plt.grid();
plt.subplot(313);
plt.plot(np.linspace(1487.03125,1322.9875-(1334.6875-1322.96875),32),np.mean(dataint,axis=1));
plt.xlabel('frequency [MHz]');
plt.ylabel('power');
plt.grid();
plt.show();
