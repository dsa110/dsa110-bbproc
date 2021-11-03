import os
import requests
import time
import pycurl
from io import BytesIO
import datetime
import math
import numpy as np
import ephem
from astropy.time import Time
import matplotlib.pyplot as plt
import dsautils.dsa_store as ds
from PyAstronomy import pyasl
import matplotlib as mpl
import matplotlib.dates as mdates
import numpy as np, matplotlib.pyplot as plt, os, pylab, glob
import scipy.signal
from scipy.stats.stats import pearsonr
from scipy import stats
from astropy.time import Time
from slack_sdk import WebClient
import os
import requests
import time
import pycurl
from io import BytesIO
import math
import ephem
import sys
sys.path.insert(0,'/mnt/data/dsa110/T3/gregtest/rfianalysis/snap_plots/')
from sat_analysis import *
from rfi_plotting import *
from opensky_api import OpenSkyApi
from mpl_toolkits.basemap import Basemap
from IPython import display
from datetime import datetime, timedelta
import mpu
import pytz
import csv
import dsautils.dsa_store as ds
mpl.rcParams['timezone'] = 'US/Pacific';


dec = np.arange(-50,120);   # declinations
ndec = len(dec);


## Time array
pat = '/home/user/vikram/PLOTS/DATA/';
dirs = glob.glob(pat+'59*');
dirs.sort();
dirs = dirs[-96:]; # = 24 hours
tim = np.linspace(float(dirs[0].split('/')[-1]),float(dirs[-1].split('/')[-1]),len(dirs)*100);


obs = ephem.Observer();
# DSA : https://github.com/dsa110/dsa110-antpos/blob/master/antpos/data/DSA110_positions_RevF.csv
obs.lat = '37.233386982';
obs.long = '-118.283405115';


syssat = ['GLONASS','GPS','GALILEO','BEIDOU'];


def gettle(system):
    file=open("/home/user/T3_detect/elev_analysis/tle.txt", 'w');
    file.write("# Last updated: %s\n"%(time.time()));

    b_obj = BytesIO();
    crl = pycurl.Curl();

    # Set URL value
    if system == 'GLONASS':
        crl.setopt(crl.URL, 'http://celestrak.com/NORAD/elements/glo-ops.txt');
    elif system == 'GPS':
        crl.setopt(crl.URL, 'http://celestrak.com/NORAD/elements/gps-ops.txt');
    elif system == 'GALILEO':
        crl.setopt(crl.URL, 'http://celestrak.com/NORAD/elements/galileo.txt');
    elif system == 'BEIDOU':
        crl.setopt(crl.URL, 'https://www.celestrak.com/NORAD/elements/beidou.txt');
    else:
        crl.setopt(crl.URL, 'https://www.celestrak.com/NORAD/elements/gnss.txt');
    crl.setopt(crl.WRITEDATA, b_obj);
    crl.perform();
    crl.close();

    # Get the content stored in the BytesIO object (in byte characters) 
    get_body = b_obj.getvalue();
    file.write(get_body.decode('utf8'));
    file.close();

ha = 0
lat = 37.233386982;

azaltall = np.zeros((2,ndec));

for s in range(len(syssat)):
    gettle(syssat[s]);
    file = open('/home/user/T3_detect/elev_analysis/tle.txt');
    tlefile = file.readlines();
    file.close();
    nSats = int((len(tlefile)-1)/3.);

    iss = [];
    for k in range(nSats):
        iss.append(ephem.readtle(tlefile[k*3+1], tlefile[k*3+2], tlefile[k*3+3]));

    possat = np.zeros((2,len(tim),nSats));
    for k in range(len(tim)):
        start_time_utc = Time(tim[k],format='mjd').to_value('datetime');
        obs.date = ephem.Date(start_time_utc);
        for kk in range(nSats):
            iss[kk].compute(obs);
            possat[0,k,kk] = (iss[kk].az);
            possat[1,k,kk] = (iss[kk].alt);

    distsat = np.zeros((ndec,len(tim),nSats))

    for k in range(ndec):
        azalt = pyasl.hadec2altaz(ha, dec[k], lat);
        obs_az = azalt[1];
        obs_el = azalt[0];
        azaltall[0,k]=obs_az;azaltall[1,k]=obs_el;
        print(syssat[s]+' -- '+str(k+1)+'/'+str(ndec));

        # if obs_el > 90.:
            # obs_el = 90.-(obs_el-90.);
            # obs_az = obs_az - 180.;

        for kk in range(nSats):
            for kkk in range(len(tim)):
                distsat[k,kkk,kk] = math.sqrt( (90.-math.degrees(possat[1,kkk,kk]))**2 + (90.-obs_el)**2 -2*(90.-math.degrees(possat[1,kkk,kk]))*(90.-obs_el)*np.cos(possat[0,kkk,kk]-math.radians(obs_az)));
    np.save('/home/user/T3_detect/elev_analysis/'+syssat[s], distsat);

np.save('/home/user/T3_detect/elev_analysis/azaltall',azaltall);


numsats5 = np.zeros((4,ndec));
numsats2 = np.zeros((4,ndec));
numsats1 = np.zeros((4,ndec));
for s in range(len(syssat)):
    distsat = np.load('/home/user/T3_detect/elev_analysis/'+syssat[s]+'.npy');
    for k in range(ndec):
        numsats5[s,k] = np.sum(distsat[k,:,:]<5.);
        numsats2[s,k] = np.sum(distsat[k,:,:]<2.);
        numsats1[s,k] = np.sum(distsat[k,:,:]<1.);

plt.figure(figsize=(10,6));
plt.plot(dec,np.sum(numsats5,axis=0)/distsat.shape[1],label='within 5 deg');
plt.plot(dec,np.sum(numsats2,axis=0)/distsat.shape[1],label='within 2 deg');
plt.plot(dec,np.sum(numsats1,axis=0)/distsat.shape[1],label='within 1 deg');
plt.legend();
plt.grid();
plt.xlabel('dec');
plt.ylabel('Prob satellite encounter');
plt.show();
