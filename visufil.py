# plot beamformed + dedispersed data
# comparison 5 beams (labeling issue)
# 2/10/21
# ghellbourg@astro.caltech.edu

import os
import os.path
import sys
import glob
import numpy as np
import yaml
import astropy.units as u
from astropy.io import ascii
import matplotlib.pyplot as plt
import json
import subprocess
from sigpyproc.Readers import FilReader

dirname  = sys.argv[1];
ibox     = int(sys.argv[2]);
specnum  = int(sys.argv[3]);
timehr   = float(sys.argv[4]);
snr      = float(sys.argv[5]);
dm       = float(sys.argv[6]);
nBeamNum = int(sys.argv[7]);

nWin = int(ibox*16);
plt.figure();
for k in range(5):
	fname = dirname + '_' + str(specnum) + '_beam' + str(nBeamNum-2+k).zfill(3) + '.fil';
	f = FilReader(fname);
	ts = f.dedisperse(dm);
	fts = ts.applyBoxcar(nWin);
	plt.plot(np.linspace(0,f.header['tsamp']*f.header['nsamples'],len(fts[nWin:-nWin])),fts[nWin:-nWin],label=str(nBeamNum-2+k));

plt.grid();
plt.legend();
plt.xlabel('time [s]');
plt.ylabel('power [arb.]');
plt.title(dirname + ' - specnum = ' + str(specnum) + ' - dedispersion = ' + str(dm) + ' - boxcarwin : ' + str(nWin));
plt.show();
