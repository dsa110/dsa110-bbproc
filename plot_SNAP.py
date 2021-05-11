# plots the latest spectra from the SNAP boards (limited to 10 SNAPs for now)

import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.time import Time

pat = '/home/user/vikram/PLOTS/DATA/';
nsnap = 10;
dirs = os.listdir(pat);
dirs.sort();
t0 = Time(float(dirs[-1]),format='mjd').to_value('datetime');
mask = np.zeros((30,2,1024));
alldata = np.zeros((30,2,1024));

x = 1530.-np.arange(1024)*250./1024.;

data = [];
ants = [];
for snap in np.arange(1,nsnap+1):

    fl = pat + dirs[-1] + '/data/snap'+str(snap)+'.npz';
    d = np.load(fl,fix_imports=True)

    ants.append(str(d['a1']))
    ants.append(str(d['a2']))
    ants.append(str(d['a3']))

    for i in range(6):
        data.append(d['specs'][i])

for i in range(len(ants)):
    y1 = np.abs(data[2*i]);
    y2 = np.abs(data[2*i+1]);
    alldata[i,0,:] = y1;
    alldata[i,1,:] = y2;

fig = plt.figure(figsize=(8,11))
plt.subplots_adjust(wspace=0.35,hspace=0.45)
for i in range(nsnap*3):
    plt.subplot(nsnap,3,i+1)
    plt.title('Antenna '+ants[i],fontsize=7.)
    if i>nsnap*3-4:
        plt.xlabel('Frequency (MHz)',fontsize=7.)
    else:
        plt.xticks([])
    plt.ylabel('power [linear]',fontsize=6.)
    plt.xticks(fontsize=7.)
    plt.yticks(fontsize=4.)

    plt.xlim(1280.,1530.)
    plt.plot(x,alldata[i,0,:])
    plt.plot(x,alldata[i,1,:],color='black',alpha=0.6);
    plt.grid();
plt.suptitle(t0.strftime('%y-%m-%d %H:%M:%S') + ' UTC');
plt.show();
