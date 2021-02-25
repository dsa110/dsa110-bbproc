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
import slack


fname = sys.argv[1];

dirname = fname[:fname.find('_')];
specnum  = int(fname[fname[fname.find('_')+1:].find('_'):][:fname[fname[fname.find('_')+1:].find('_'):].find('_')]);

jsonfile = glob.glob('/mnt/data/dsa110/T3/corr00/' + dirname + '/*'+str(specnum)+'.json')[0];
with open(jsonfile) as f:
    data = json.load(f);
timehr   = float(data.get(list(data.keys())[0]).get('mjds'));
snr      = float(data.get(list(data.keys())[0]).get('snr'));
dm       = float(data.get(list(data.keys())[0]).get('dm'));
nBeamNum = int(data.get(list(data.keys())[0]).get('ibeam'));
ibox     = int(data.get(list(data.keys())[0]).get('ibox'));


nWin = int(ibox*16);
f = FilReader(fname);
ts = f.dedisperse(dm);
fts = ts.applyBoxcar(nWin);

plt.figure();
plt.plot(np.linspace(0,f.header['tsamp']*f.header['nsamples'],len(fts[nWin:-nWin])),fts[nWin:-nWin]);
plt.grid();
plt.xlabel('time [s]');
plt.ylabel('power [arb.]');
plt.title(dirname + ' - specnum = ' + str(specnum) + ' - dedispersion = ' + str(dm) + ' - boxcarwin : ' + str(nWin));
plt.savefig(fname.replace('.fil','.png'));

client = slack.WebClient(token='xoxb-508911196752-1791739181778-i0912FjZnPYtLJnKIlyKSeHN');
client.files_upload(channels='candidates',file=fname.replace('.fil','.png'),initial_comment=fname.replace('.fil','.png'))
#plt.show();
