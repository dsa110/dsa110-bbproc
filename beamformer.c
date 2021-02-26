# plot triggered FRB candidates
# liam.dean.connor@gmail.com & ghellbourg@astro.caltech.edu
# 25/02/2021

import os
import os.path
import sys


import numpy as np
import matplotlib.pyplot as plt 
import json
import glob

#import filterbank
from sigpyproc.Readers import FilReader
import slack

plt.rcParams.update({
                    'font.size': 12,
                    'font.family': 'serif',
                    'axes.labelsize': 14,
                    'axes.titlesize': 15,
                    'xtick.labelsize': 12,
                    'ytick.labelsize': 12,
                    'xtick.direction': 'in',
                    'ytick.direction': 'in',
                    'xtick.top': True,
                    'ytick.right': True,
                    'lines.linewidth': 0.5,
                    'lines.markersize': 5,
                    'legend.fontsize': 14,
                    'legend.borderaxespad': 0,
                    'legend.frameon': False,
                    'legend.loc': 'lower right'})


def read_fil_data_dsa(fn, start=0, stop=1):
    """ Read in filterbank data
    """
    fil_obj = FilReader(fn)
    header = fil_obj.header
    delta_t = fil_obj.header['tsamp'] # delta_t in seconds                                                                                                                  
    fch1 = header['fch1']
    nchans = header['nchans']
    foff = header['foff']
    fch_f = fch1 + nchans*foff
    freq = np.linspace(fch1,fch_f,nchans)

    try:
        data = fil_obj.readBlock(start, stop)
    except(ValueError):
        data = 0

    return data, freq, delta_t, header

def plotfour(dataft, datats, datadmt, 
             multibeam=None, figname_out=None, dm=0,
             dms=[0,1], 
             datadm0=None, suptitle='', heimsnr=-1,
             ibox=1, ibeam=-1, prob=-1):
    """ Plot a trigger's dynamics spectrum, 
        dm/time array, pulse profile, 
        multibeam info (optional), and zerodm (optional)

        Parameter
        ---------
        dataft : 
            freq/time array (nfreq, ntime)
        datats : 
            dedispersed timestream
        datadmt : 
            dm/time array (ndm, ntime)
        multibeam : 
            TBD
        figname_out : 
            save figure with this file name 
        dm : 
            dispersion measure of trigger 
        dms : 
            min and max dm for dm/time array 
        datadm0 : 
            raw data timestream without dedispersion
    """
    datats /= np.std(datats[datats!=np.max(datats)])
    nfreq, ntime = dataft.shape
    dm_min, dm_max = dms[0], dms[1]
    tmin, tmax = 1e3*dataft.header['tsamp'], 1e3*dataft.header['tsamp']*ntime
#    freqmin, freqmax = min(freq), max(freq)
    freqmax = dataft.header['fch1']
    freqmin = freqmax + dataft.header['nchans']*dataft.header['foff']
    tarr = np.linspace(tmin, tmax, ntime)
    fig = plt.figure(figsize=(8,10))

    plt.subplot(321)
    extentft=[tmin,tmax,freqmax,freqmin]
    plt.imshow(dataft, aspect='auto',extent=extentft)
    plt.xlabel('Time (ms)')
    plt.ylabel('Freq (MHz)')

    plt.subplot(322)
    extentdm=[tmin, tmax, dm_min, dm_max]
    plt.imshow(datadmt, aspect='auto',extent=extentdm)
    plt.xlabel('Time (ms)')
    plt.ylabel(r'DM (pc cm$^{-3}$)')

    plt.subplot(323)
    plt.plot(tarr, datats)
    plt.grid('on', alpha=0.25)
    plt.xlabel('Time (ms)')
    plt.ylabel(r'Power ($\sigma$)')
    plt.text(0.6*(tmin+tmax), 0.5*(max(datats)+np.median(datats)), 
            'Heimdall S/N : %0.1f\nHeimdall DM : \
            %d\nHeimdall ibox : %d\nibeam : %d\nprob : %0.2f' % (heimsnr,dm,ibox,ibeam,prob), 
            fontsize=8, verticalalignment='center')
    plt.subplot(324)
    plt.xticks([])
    plt.yticks([])

    if multibeam is None:
        plt.text(0.20, 0.55, 'Multibeam info\nunder construction',
                fontweight='bold')

    # if figname_out is not None:
        # plt.savefig(figname_out)

    if datadm0 is not None:
        # tobs = datadm0.header['tobs']
        # tsdm0 = datadm0.mean(0)
        plt.subplot(313)
#        plt.plot(np.linspace(0, tobs, len(tsdm0)), tsdm0, color='k')
        plt.plot(np.linspace(0, tmax, len(datadm0)), datadm0, color='k')
        plt.legend(['DM=0 Timestream'], loc=2)
        plt.xlabel('Time (ms)')

    plt.suptitle(suptitle, c='C1')
    plt.tight_layout()
    #plt.show()
    if figname_out is not None:
        plt.savefig(figname_out)

def dm_transform(data, dm_max=20,
                 dm_min=0, dm0=None, ndm=64, 
                 freq_ref=None, downsample=16):
    """ Transform freq/time data to dm/time data.                                                                                                                                           
    """
    ntime = data.shape[1]

    dms = np.linspace(dm_min, dm_max, ndm, endpoint=True)

    if dm0 is not None:
        dm_max_jj = np.argmin(abs(dms-dm0))
        dms += (dm0-dms[dm_max_jj])

    data_full = np.zeros([ndm, ntime//downsample])

    for ii, dm in enumerate(dms):
        dd = data.dedisperse(dm)
#        print(dd.dm)
        _dts = np.mean(dd,axis=0)
        data_full[ii] = _dts[:ntime//downsample*downsample].reshape(ntime//downsample, downsample).mean(1)

    return data_full, dms

def proc_cand_fil(fnfil, dm, ibox, snrheim=-1, 
                  pre_rebin=8, nfreq_plot=64,
                  heim_raw_tres=32, 
                  rficlean=False, ndm=64):
    """ Take filterbank file path, preprocess, and 
    plot trigger

    Parameters:
    ----------

    fnfil   : str 
        path to .fil file 
    DM      : float 
        dispersion measure of trigger 
    ibox    : int 
        preferred boxcar width 
    snrheim : float 
        S/N of candidate found by Heimdall
    pre_rebin : int 
        rebin in time by this factor *before* dedispersion (saves time)
    nfreq_plot : int 
        number of frequency channels in output
    heim_raw_tres : 32  
    """
    data, freq, delta_t_raw, header = read_fil_data_dsa(fnfil, start=0, 
                                                       stop=int(150000))
    nfreq0, ntime0 = data.shape

    # Rebin in frequency by 4x
    data = data.reshape(nfreq0//4, 4, ntime0).mean(1)
    data.header['foff'] = data.header['foff']*4
    data = data.downsample(pre_rebin)

    if rficlean:
#        print("Cleaning data perchannel")
        data = cleandata(data, clean_type='perchannel')

    tsdm0 = np.mean(data,axis=0)
#    data = data.dedisperse(dm)

    datadm, dms = dm_transform(data, dm_max=dm+250,
                               dm_min=dm-250, dm0=dm, ndm=ndm, 
                               freq_ref=None, 
                               downsample=heim_raw_tres*ibox)#//pre_rebin)
    data = data.dedisperse(dm)
    data = data.downsample(heim_raw_tres*ibox)#//pre_rebin)
    data = data.reshape(nfreq_plot, data.shape[0]//nfreq_plot, 
                        data.shape[1]).mean(1)

    data = data-np.median(data,axis=1,keepdims=True)
    data /= np.std(data)

    return data, datadm, tsdm0, dms


def cleandata(data, threshold_time=3.25, threshold_frequency=2.75, bin_size=32,
              n_iter_time=3, n_iter_frequency=3, clean_type='time', wideclean=None):
    """ Take filterbank object and mask
    RFI time samples with average spectrum.
    Parameters:
    ----------
    data :
        data array (nfreq, ntime)
    threshold_time : float
        units of sigma
    threshold_frequency : float
        units of sigma
    bin_size : int
        quantization bin size
    n_iter_time : int
        Number of iteration for time cleaning
    n_iter_frequency : int
        Number of iteration for frequency cleaning
    clean_type : str
        type of cleaning to be done.
        Accepted values: 'time', 'frequency', 'both', 'perchannel'
    Returns:
    -------
    cleaned filterbank object
    """
    if clean_type not in ['time', 'both', 'frequency', 'perchannel']:
        return data
        
    nfreq = data.shape[0]

    dtmean = np.mean(data, axis=-1)
    # Clean in time
    #sys_temperature_bandpass(data.data)
    #remove_noisy_freq(data.data, 3)
    #remove_noisy_channels(data.data, sigma_threshold=2, iters=5)
    if clean_type in ['time', 'both']:
        for i in range(n_iter_time):
            dfmean = np.mean(data, axis=0)
            stdevf = np.std(dfmean)
            medf = np.median(dfmean)
            maskf = np.where(np.abs(dfmean - medf) > threshold_time*stdevf)[0]
            # replace with mean spectrum
            data[:, maskf] = dtmean[:, None]*np.ones(len(maskf))[None]

    if clean_type=='perchannel':
        for ii in range(n_iter_time):
            dtmean = np.mean(data, axis=1, keepdims=True)
            dtsig = np.std(data, axis=1)
            for nu in range(data.shape[0]):
                d = dtmean[nu]
                sig = dtsig[nu]
                maskpc = np.where(np.abs(data[nu]-d)>threshold_time*sig)[0]
                data[nu][maskpc] = d

    # Clean in frequency
    # remove bandpass by averaging over bin_size ajdacent channels
    if clean_type in ['frequency', 'both']:
        for ii in range(n_iter_frequency):
            dtmean_nobandpass = data.mean(1) - dtmean.reshape(-1, bin_size).mean(-1).repeat(bin_size)
            stdevt = np.std(dtmean_nobandpass)
            medt = np.median(dtmean_nobandpass)
            maskt = np.abs(dtmean_nobandpass - medt) > threshold_frequency*stdevt
            data[maskt] = np.median(dtmean)#dtmean.reshape(-1, bin_size).mean(-1).repeat(bin_size)[maskt]

    return data

def plot_fil(fn, dm, ibox, multibeam=None, figname_out=None,
             ndm=32, suptitle='', heimsnr=-1,
             ibeam=-1, rficlean=True, nfreq_plot=64, 
             classify=False):

    dataft, datadm, tsdm0, dms = proc_cand_fil(fn, dm, ibox, snrheim=-1, 
                  pre_rebin=8, nfreq_plot=nfreq_plot, ndm=ndm, rficlean=rficlean,
                  heim_raw_tres=32)

    if classify:
        pass
    else:
        prob = -1

    plotfour(dataft, dataft.mean(0), datadm, datadm0=tsdm0, 
             multibeam=None, figname_out=figname_out, dm=dm,
             dms=[dms[0],dms[-1]], 
             suptitle=suptitle, heimsnr=heimsnr,
             ibox=ibox, ibeam=ibeam, prob=prob)


if __name__=='__main__':

    fname = sys.argv[1];
    dirname = fname[:fname.find('_')];
    specnum  = int(fname[fname[fname.find('_')+1:].find('_'):][:fname[fname[fname.find('_')+1:].find('_'):].find('_')])

    jsonfile = glob.glob('/mnt/data/dsa110/T3/corr00/' + dirname + '/*'+str(specnum)+'.json')[0]

    with open(jsonfile) as f:
        triggerdata = json.load(f)

    timehr   = float(triggerdata.get(list(triggerdata.keys())[0]).get('mjds'))
    snr      = float(triggerdata.get(list(triggerdata.keys())[0]).get('snr'))
    dm       = float(triggerdata.get(list(triggerdata.keys())[0]).get('dm'))
    nBeamNum = int(triggerdata.get(list(triggerdata.keys())[0]).get('ibeam'))
    ibox     = int(triggerdata.get(list(triggerdata.keys())[0]).get('ibox'))
    nWin = ibox*32;

    suptitle = dirname + ' - specnum = ' + str(specnum) \
              + ' - dedispersion = ' + str(dm) + ' - boxcarwin : ' + str(nWin);

    fnameout = fname.replace('.fil','.png');
    plot_fil(fname, dm, ibox, multibeam=None, figname_out=fnameout,
                 ndm=32, suptitle=suptitle, heimsnr=snr,
                 ibeam=nBeamNum, rficlean=True, nfreq_plot=64, 
                 classify=False);
    client = slack.WebClient(token='**********************');
    client.files_upload(channels='candidates',file=fnameout,initial_comment=fnameout);
