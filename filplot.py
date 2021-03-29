# plot triggered FRB candidates
# liam.dean.connor@gmail.com & ghellbourg@astro.caltech.edu
# 25/02/2021

import os
import os.path
import sys

import scipy.signal
from scipy import stats

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import json
import glob
import optparse

#import filterbank
from sigpyproc.Readers import FilReader
import slack

# Keras neural network model for Freq/Time array
MLMODELPATH='/home/user/connor/software/machine_learning/20190501freq_time.hdf5'
BASEDIR='/mnt/data/dsa110/'

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
             beam_time_arr=None, figname_out=None, dm=0,
             dms=[0,1],
             datadm0=None, suptitle='', heimsnr=-1,
             ibox=1, ibeam=-1, prob=-1, showplot=True):
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
        beam_time_arr :
            beam time SNR array (nbeam, ntime)
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
    tmin, tmax = 0, 1e3*dataft.header['tsamp']*ntime
    freqmax = dataft.header['fch1']
    freqmin = freqmax + dataft.header['nchans']*dataft.header['foff']
    tarr = np.linspace(tmin, tmax, ntime)
    fig = plt.figure(figsize=(8,10))

    plt.subplot(321)
    extentft=[tmin,tmax,freqmin,freqmax]
    plt.imshow(dataft, aspect='auto',extent=extentft, interpolation='nearest')
    plt.xlim(0,1000)
    plt.xlabel('Time (ms)')
    plt.ylabel('Freq (MHz)')
    if prob!=-1:
        plt.text(tmin+50,0.5*(freqmax+freqmin),"Prob=%0.2f" % prob, color='white', fontweight='bold')
    plt.subplot(322)
    extentdm=[tmin, tmax, dm_max, dm_min]
    plt.imshow(datadmt, aspect='auto',extent=extentdm)
    plt.xlim(0,1000)
    plt.xlabel('Time (ms)')
    plt.ylabel(r'DM (pc cm$^{-3}$)')

    plt.subplot(323)
    plt.plot(tarr, datats)
    plt.grid('on', alpha=0.25)
    plt.xlabel('Time (ms)')
    plt.ylabel(r'Power ($\sigma$)')
    plt.xlim(0,1000)
    plt.text(0.55*(tmin+1000.), 0.5*(max(datats)+np.median(datats)),
            'Heimdall S/N : %0.1f\nHeimdall DM : %d\
            \nHeimdall ibox : %d\nibeam : %d' % (heimsnr,dm,ibox,ibeam),
            fontsize=8, verticalalignment='center')
    plt.subplot(324)


    if beam_time_arr is None:
        plt.xticks([])
        plt.yticks([])
        plt.text(0.20, 0.55, 'Multibeam info\nunder construction',
                fontweight='bold')
    else:
        plt.imshow(beam_time_arr, aspect='auto', extent=[tmin, tmax, beam_time_arr.shape[0], 0],
                  interpolation='nearest')
        plt.axvline(540, ymin=0, ymax=6, color='r', linestyle='--', alpha=0.55)
        plt.axvline(460, ymin=0, ymax=6, color='r', linestyle='--', alpha=0.55)
        plt.axhline(max(0,ibeam-1), xmin=0, xmax=100, color='r', linestyle='--', alpha=0.55)
        plt.axhline(ibeam+1, xmin=0, xmax=100, color='r', linestyle='--', alpha=0.55)
        plt.xlim(0, 1000)
        plt.xlabel('Time (ms)')
        plt.ylabel('Beam', fontsize=15)

    # if figname_out is not None:
        # plt.savefig(figname_out)

    if datadm0 is not None:
        # tsdm0 = datadm0.mean(0)
        plt.subplot(325)
#        plt.plot(np.linspace(0, tobs, len(tsdm0)), tsdm0, color='k')
        plt.plot(np.linspace(0, tmax, len(datadm0[0])), datadm0.mean(0), color='k')
        plt.legend(['DM=0 Timestream'], loc=2, fontsize=10)
        plt.xlabel('Time (ms)')

        plt.subplot(326)
        plt.plot(np.linspace(freqmax,freqmin,datadm0.shape[0]), np.mean(datadm0,axis=-1), color='k')
        plt.semilogy()
        plt.legend(['spectrum'], loc=2)
        plt.xlabel('freq [MHz]')

    plt.suptitle(suptitle, c='C1')
    plt.tight_layout()
    if figname_out is not None:
        plt.savefig(figname_out)
    if showplot:
        plt.show()


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
        _dts = np.mean(dd,axis=0)
        data_full[ii] = _dts[:ntime//downsample*downsample].reshape(ntime//downsample, downsample).mean(1)

    return data_full, dms

def proc_cand_fil(fnfil, dm, ibox, snrheim=-1,
                  pre_rebin=1, nfreq_plot=64,
                  heim_raw_tres=1,
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
    header = read_fil_data_dsa(fnfil, 0, 1)[-1]
    # read in 4 seconds of data
    nsamp = int(4.0/header['tsamp'])
    data, freq, delta_t_raw, header = read_fil_data_dsa(fnfil, start=0,
                                                       stop=nsamp)
    nfreq0, ntime0 = data.shape

    # Ensure that you do not pre-downsample by more than the total boxcar
    pre_rebin = min(pre_rebin, ibox*heim_raw_tres)

    # Rebin in frequency by 8x
    data = data.reshape(nfreq0//8, 8, ntime0).mean(1)   #bin 8 channels
    data.header['foff'] = data.header['foff']*8         #bin 8 channels
    data.header['nchans'] = data.header['nchans']/8
    data = data.downsample(pre_rebin)

    datadm0 = data.copy()

    if rficlean:
#        print("Cleaning data perchannel")
        data = cleandata(data, clean_type='aladsa')

    tsdm0 = np.mean(data,axis=0)

    datadm, dms = dm_transform(data, dm_max=dm+250,
                               dm_min=dm-250, dm0=dm, ndm=ndm,
                               freq_ref=None,
                               downsample=heim_raw_tres*ibox//pre_rebin)
    data = data.dedisperse(dm)
    data = data.downsample(heim_raw_tres*ibox//pre_rebin)
    data = data.reshape(nfreq_plot, data.shape[0]//nfreq_plot,
                        data.shape[1]).mean(1)

    data = data-np.median(data,axis=1,keepdims=True)
    data /= np.std(data)

    return data, datadm, tsdm0, dms, datadm0


def medflagdata(spec, filtsize, thres):
    specfilt = scipy.signal.medfilt(spec,kernel_size=int(filtsize));
    speccorrec = spec - specfilt;
    specstd = stats.median_absolute_deviation(speccorrec);
    return np.concatenate((np.argwhere(speccorrec > thres*specstd),np.argwhere(speccorrec < -thres*specstd)))

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
    if clean_type not in ['time', 'both', 'frequency', 'perchannel', 'aladsa']:
        return data

    nfreq = data.shape[0]
    ntimes = data.shape[1]

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

    if clean_type=='aladsa':
        print('flagging a la DSA\n');
        meanidx = medflagdata(dtmean, 21, 5.);
        varidx = medflagdata(np.var(data,axis=-1), 21, 5.);
        allidx = np.concatenate((meanidx,varidx));
        allidx = np.asarray(list(set(list(np.ravel(allidx)))));
        data[allidx,:] = np.zeros((len(allidx),ntimes));


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

def generate_beam_time_arr(fl, ibeam=0, pre_rebin=1,
                           dm=0, ibox=1, heim_raw_tres=1):
    """ Take list of nbeam .fil files, dedisperse each
    to the dm of the main trigger, and generate an
    (nbeam, ntime) SNR array.

    Parameters:
    -----------
    fl : list
        list of .fil files, each 4 seconds long
    ibeam : int
        beam number of trigger
    pre_rebin :
        downsample by this factor before dedispersion to save time
    dm : int
        dm of ibeam candidate
    ibox : int
        boxcar width of ibeam candidate
    heim_raw_tres : int
        ratio of

    Returns:
    --------
    beam_time_arr : ndarray
        array of SNR values (nbeam, ntime)
    """
    fl.sort()
    nbeam = len(fl)
    header = read_fil_data_dsa(fl[0], 0, 1)[-1]
    # read in 4 seconds of data
    nsamp = int(4.0/header['tsamp'])
    nsamp_final = nsamp // (heim_raw_tres*ibox)

    beam_time_arr = np.zeros([nbeam, nsamp_final])

    for fnfil in fl:
        print(fnfil, beam_time_arr.shape)
        beamno = int(fnfil.strip('.fil').split('_')[-1])
        data, freq, delta_t_raw, header = read_fil_data_dsa(fnfil, start=0,
                                                           stop=nsamp)
        nfreq0, ntime0 = data.shape

        # Ensure that you do not pre-downsample by more than the total boxcar
        pre_rebin = min(pre_rebin, ibox*heim_raw_tres)

        # Rebin in frequency by 8x
        data = data.reshape(nfreq0//8, 8, ntime0).mean(1)   #bin 8 channels
        data.header['foff'] = data.header['foff']*8         #bin 8 channels
        data.header['nchans'] = data.header['nchans']/8
        data = data.downsample(pre_rebin)
        data = data.dedisperse(dm)
        data = data.downsample(heim_raw_tres*ibox//pre_rebin)
        datats = np.mean(data, axis=0)

        # Normalize data excluding outliers
        datatscopy = datats.copy()
        datatscopy.sort()
        medts = np.median(datatscopy[:int(0.975*len(datatscopy))])
        sigts = np.std(datatscopy[:int(0.975*len(datatscopy))])
        datats -= medts
        datats /= sigts

        beam_time_arr[beamno, :] = datats

    return beam_time_arr


def plot_fil(fn, dm, ibox, multibeam=None, figname_out=None,
             ndm=32, suptitle='', heimsnr=-1,
             ibeam=-1, rficlean=True, nfreq_plot=32,
             classify=False, heim_raw_tres=1, showplot=True):
    """ Vizualize FRB candidates on DSA-110
    """
    if type(multibeam)==list:
        beam_time_arr = generate_beam_time_arr(multibeam, ibeam=ibeam, pre_rebin=1,
                                             dm=dm, ibox=ibox,
                                             heim_raw_tres=heim_raw_tres)
    else:
        beam_time_arr = None


    dataft, datadm, tsdm0, dms, datadm0 = proc_cand_fil(fn, dm, ibox, snrheim=-1,
                                               pre_rebin=1, nfreq_plot=nfreq_plot,
                                               ndm=ndm, rficlean=rficlean,
                                               heim_raw_tres=heim_raw_tres)

    if classify:
        from keras.models import load_model
        fnmodel=MLMODELPATH
        model = load_model(fnmodel)
        mm = np.argmax(dataft.mean(0))
        tlow, thigh = mm-32, mm+32
        if mm<32:
            tlow=0
            thigh=64
        if thigh>dataft.shape[1]:
            thigh=dataft.shape[1]
            tlow=thigh-64
        dataml = dataft[:,tlow:thigh]
        dataml -= np.median(dataml, axis=1, keepdims=True)
        dataml /= np.std(dataml, axis=-1)[:, None]
        dataml[dataml!=dataml] = 0.0
        dataml = dataml[None,..., None]
        prob = model.predict(dataml)[0,1]
    else:
        prob = -1

    plotfour(dataft, dataft.mean(0), datadm, datadm0=datadm0,
             beam_time_arr=beam_time_arr, figname_out=figname_out, dm=dm,
             dms=[dms[0],dms[-1]],
             suptitle=suptitle, heimsnr=heimsnr,
             ibox=ibox, ibeam=ibeam, prob=prob, showplot=showplot)

def read_json(jsonfile):
    with open(jsonfile) as f:
        triggerdata = json.load(f)

    timehr   = float(triggerdata.get(list(triggerdata.keys())[0]).get('mjds'))
    snr      = float(triggerdata.get(list(triggerdata.keys())[0]).get('snr'))
    dm       = float(triggerdata.get(list(triggerdata.keys())[0]).get('dm'))
    ibeam = int(triggerdata.get(list(triggerdata.keys())[0]).get('ibeam'))
    ibox     = int(triggerdata.get(list(triggerdata.keys())[0]).get('ibox'))

    return timehr,snr,dm,ibeam,ibox


if __name__=='__main__':

    parser = optparse.OptionParser(prog="filplotter",
                                   version="",
                                   usage="%prog fname [OPTIONS]",
                                   description="Visualize and classify filterbank data")

    parser.add_option('-s', '--slack', dest='slack', action="store_true",help="send figure to slack")
    parser.add_option('-d', '--dm', dest='dm',
                      help="DM ", default=None)
    parser.add_option('-c', '--classify', dest='classify', action="store_true",
                      help="classify using ML")
    parser.add_option('-r', '--rficlean', dest='rficlean', action="store_true",
                      help="excise RFI from data")
    parser.add_option('-w', '--ibox', dest='ibox', type=int,
                      help="ibox found by Heimdall", default=1)
    parser.add_option('--ndm', dest='ndm', type=int, default=32,
                      help="number of DMs for DM/time plot")
    parser.add_option('--ntime_plot', dest='ntime_plot', type=int, default=64,
                      help="number of samples to plot")
    parser.add_option('--nfreq_plot', dest='nfreq_plot', type=int, default=32,
                      help="number of freq channels to plot")

    options, args = parser.parse_args()
    datestr = args[0]
    specnum = args[1]

    flist = glob.glob(BASEDIR+'/T1/corr*/'+datestr+'/fil_%s/*.fil' % specnum)
    jsonfile = glob.glob(BASEDIR+'/T3/corr01/'+datestr+'/*%s*.json' % specnum)[0]

    timehr,snr,dm,ibeam,ibox = read_json(jsonfile)
    print('Read JSON file')

    beamindlist=[]
    for fnfil in flist:
        beamno = int(fnfil.strip('.fil').split('_')[-1])
        beamindlist.append(beamno)
        if beamno==ibeam:
            fname = fnfil

    if options.slack:
        showplot=False
    else:
        showplot=True

    if options.dm is not None:
        dm = float(options.dm)
        ibox = int(options.ibox)
        plot_fil(fname, dm, ibox, multibeam=None, figname_out=None,
                 ndm=options.ndm, suptitle='', heimsnr=-1,
                 ibeam=1, rficlean=options.rficlean, nfreq_plot=options.nfreq_plot,
                 classify=options.classify, heim_raw_tres=1, showplot=showplot)
        exit()

    # fi = fname.find('_');
    # se = fname[fi+1:].find('_');
    # dirname = fname[:fi];
    # specnum  = int(fname[fi+1:fi+se+1])
    # print('dirname = '+dirname+' -- specnum = '+str(specnum));

#    dirname = fname[:fname.find('_')];
#    specnum  = int(fname[fname[fname.find('_')+1:].find('_'):][:fname[fname[fname.find('_')+1:].find('_'):].find('_')])

    # For now this will be a placeholder.
#    jsonfile = glob.glob('/mnt/data/dsa110/T3/corr00/' + dirname + '/*'+str(specnum)+'.json')[0]

    suptitle = 'specnum:%s  dm:%0.4f  boxcar:%d  ibeam:%d' % (specnum, dm, int(ibox), int(ibeam))

    fnameout = fname.replace('.fil','.png')

    plot_fil(fname, dm, ibox, figname_out=fnameout,
             ndm=options.ndm, suptitle=suptitle, heimsnr=snr,
             ibeam=ibeam, rficlean=options.rficlean,
             nfreq_plot=options.nfreq_plot,
             classify=options.classify, showplot=showplot,
             multibeam=flist, heim_raw_tres=1)

    if options.slack:
        print("Sending to slack")
        client = slack.WebClient(token='XXXXXXXXXXXXXXXX');
        client.files_upload(channels='candidates',file=fnameout,initial_comment=fnameout);
