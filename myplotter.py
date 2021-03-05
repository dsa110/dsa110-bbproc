# takes a FIL file as argument
# plots de-dispersed profile with / without RFI mitigation
# Also plots freq spectrum and highlights flagged out channels
# posts plot to slack

import os
import os.path
import sys
import scipy.signal
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt 
import json
import glob
from sigpyproc.Readers import FilReader
import slack

def medflagdata(spec, filtsize, thres):
    specfilt = scipy.signal.medfilt(spec,kernel_size=int(filtsize));
    speccorrec = spec - specfilt;
    specstd = stats.median_absolute_deviation(speccorrec);
    return np.concatenate((np.argwhere(speccorrec > thres*specstd),np.argwhere(speccorrec < -thres*specstd)));

fname = sys.argv[1];
#fname = '/mnt/data/dsa110/findpulse/28feb21_24230272_1514869_beam163.fil'
fn = fname.split('/')[-1];
fi = fn.find('_');
se = fn[fi+1:].find('_');
dirname = fn[:fi];
specnum  = int(fn[fi+1:fi+se+1]);

jsonfile = glob.glob('/mnt/data/dsa110/T3/corr00/' + dirname + '/*'+str(specnum)+'.json')[0]

with open(jsonfile) as f:
    triggerdata = json.load(f)

timehr   = float(triggerdata.get(list(triggerdata.keys())[0]).get('mjds'))
snr      = float(triggerdata.get(list(triggerdata.keys())[0]).get('snr'))
dm       = float(triggerdata.get(list(triggerdata.keys())[0]).get('dm'))
nBeamNum = int(triggerdata.get(list(triggerdata.keys())[0]).get('ibeam'))
ibox     = int(triggerdata.get(list(triggerdata.keys())[0]).get('ibox'))

# read data, downsample, dedisperse
fil_obj = FilReader(fname);
delta_t = fil_obj.header['tsamp']; # delta_t in seconds                                                                                                                  
data = fil_obj.readBlock(0, -1);
data = data.reshape(data.shape[0]//8, 8, data.shape[1]).mean(1)   #bin 8 channels
#data = data.downsample(32);
data /= np.std(data,axis=1,keepdims=True);
data.header['foff'] = data.header['foff']*8         #bin 8 channels
data.header['nchans'] = data.header['nchans'] / 8;
#data.header['tsamp'] = data.header['tsamp']*32.;
data = data.dedisperse(dm);

dtmean = np.mean(data,axis=-1);
meanidx = medflagdata(dtmean, 21, 5.);
varidx = medflagdata(np.var(data,axis=-1), 21, 5.);
allidx = np.concatenate((meanidx,varidx));
allidx = np.asarray(list(set(list(np.ravel(allidx)))));
#data = data - np.mean(data,axis=-1,keepdims=True);
dataflagged = np.copy(data);
dataflagged[allidx,:] = np.zeros((len(allidx),data.shape[1]));

nWin = int(np.floor(ibox*32./2.));

tstr1 = np.mean(data,axis=0);
tstr1 = np.convolve(tstr1,np.ones(nWin));

tstr2 = np.mean(dataflagged,axis=0);
tstr2 = np.convolve(tstr2,np.ones(nWin));

fch1 = data.header['fch1'];
foff = data.header['foff'];
plt.figure();
plt.suptitle(fname.split('/')[-1]+' - DM='+str(dm))
plt.subplot(121);
plt.plot((np.arange(len(tstr1[nWin:-nWin])))*data.header['tsamp'],tstr1[nWin:-nWin]);
plt.plot((np.arange(len(tstr1[nWin:-nWin])))*data.header['tsamp'],tstr2[nWin:-nWin]);
plt.xlabel('time [s]');
plt.grid();
plt.subplot(122);
plt.plot(np.linspace(fch1,fch1+foff*data.shape[0],data.shape[0]),np.mean(data,axis=-1));
plt.plot(np.linspace(fch1,fch1+foff*data.shape[0],data.shape[0]),np.mean(dataflagged,axis=-1),'.r');
plt.xlabel('frequency [MHz]');
plt.grid();
#plt.show();

pngname = '/mnt/data/dsa110/findpulse/codes/plots/'+fn.replace('.fil','.png');
plt.savefig(pngname);
client = slack.WebClient(token='XXXXXXXXXXXXX');
client.files_upload(channels='UNDTK3TLN',file=pngname,initial_comment=pngname.split('/')[-1]);
