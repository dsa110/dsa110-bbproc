import os
import tqdm
import h5py
import numpy as np
from pkg_resources import resource_filename
import astropy.units as u
import astropy.constants as c
from antpos.utils import get_itrf
from dsacalib.fringestopping import calc_uvw
import dsacalib.constants as ct
from dsautils import cnf
from astropy.time import Time
import json, sys

def get_cable_delays(params: dict) -> np.ndarray:
    """Calculate cable delays from the measured outrigger cable delays.

    Antenna names in both input parameters are indexed at 1.

    Parameters
    ----------
    params : dict
        vis params.

    Returns
    -------
    np.ndarray
        The cable delay for each baseline in bname.
    """
    delays = np.zeros(len(params['bname']), dtype=np.int)
    for i, bn in enumerate(params['bname']):
        ant1, ant2 = bn.split('-')
        delays[i] = params['antenna_cable_delays'].get(int(ant1), 0)-\
                    params['antenna_cable_delays'].get(int(ant2), 0)
    return delays

def calculate_uvw_and_geodelay(vis_params: dict, pt_dec: float) -> tuple:
    """Calculate the uvw coordinates in the correlated file.
    
    Parameters
    ----------
    vis_params : dict
        Parameters describing the correlated visibilities.
    tobs : float
        The time in mjd.
    pt_dec : float
        The pointing declination of the observation in rad.

    Returns
    -------
    buvw : np.ndarray(float)
        The uvw coordinates of each baseline, dimensions (1, baseline, 3).
    ant_bw : np.ndarray(float)
        The geometric path length, w, to each antenna, dimensions (1, antenna).
    """
    bu, bv, bw = calc_uvw(vis_params['blen'], np.asarray([vis_params['tobs']]), 'HADEC',
                          np.tile(0.*u.rad, 1), np.tile(pt_dec, 1))
    buvw = np.array([bu, bv, bw]).T

    ant_bw = buvw[:, vis_params['refidxs'], -1]

    return buvw, ant_bw


def get_blen(antennas: list) -> tuple:
    """Gets the baseline lengths for a subset of antennas.

    Parameters
    ----------
    antennas : list
        The antennas used in the array.

    Returns
    -------
    blen : array
        The ITRF coordinates of all of the baselines.
    bname : list
        The names of all of the baselines.
    """
    ant_itrf = get_itrf(
        latlon_center=(ct.OVRO_LAT*u.rad, ct.OVRO_LON*u.rad, ct.OVRO_ALT*u.m)
    ).loc[antennas]
    xx = np.array(ant_itrf['dx_m'])
    yy = np.array(ant_itrf['dy_m'])
    zz = np.array(ant_itrf['dz_m'])
    # Get uvw coordinates
    nants = len(antennas)
    nbls = (nants*(nants+1))//2
    blen = np.zeros((nbls, 3))
    bname = []
    k = 0
    for i in range(nants):
        for j in range(i, nants):
            blen[k, :] = np.array([
                xx[i]-xx[j],
                yy[i]-yy[j],
                zz[i]-zz[j]
            ])
            bname += ['{0}-{1}'.format(
                antennas[i],
                antennas[j]
            )]
            k += 1
    return blen, bname


def gen_vis_params(tstart: 'astropy.time.Time'):
    
    conf = cnf.Conf()
    corrconf = conf.get('corr')
    mfsconf = conf.get('fringe')

    antenna_order = list(corrconf['antenna_order'].values())[:63]
    outrigger_delays = mfsconf['outrigger_delays']
    nant = len(antenna_order)
    nbls = (nant*(nant+1))//2
    blen, bname = get_blen(antenna_order)
    refidxs = []
    refant = str(antenna_order[0])
    for i, bn in enumerate(bname):
        if refant in bn.split('-'):
            refidxs += [i]
    vis_params = {
        # Baselines
        'antenna_order': antenna_order,
        'blen': blen,
        'bname': bname,
        'nbls': nbls,
        'refidxs': refidxs,
        'antenna_cable_delays': outrigger_delays,
        # Time
        'tobs': tstart.mjd,
    }

    cable_delays = get_cable_delays(vis_params)
    vis_params['baseline_cable_delays'] = cable_delays
    
    return vis_params
    
    
def get_total_delay(ant_bw: np.ndarray, vis_params: dict) -> np.ndarray:
    """Calculate total (cable plus geometric) delay for each baseline.

    Parameters
    ----------
    baseline_cable_delay : np.ndarray
        The relative cable delay on each baseline to be accounted for.
    ant_bw : np.ndarray
        The w-term for each antenna, describing the geometric delay.
    vis_params : dict
        The vis params.

    Returns
    -------
    np.ndarray
        The total delay for each baseline, including cable and geometric terms.
        Dimensions (baseline).
    """
    nbls = len(vis_params['bname'])
    ntimes = ant_bw.shape[0]
    total_delay = np.zeros((ntimes, nbls))
    for bni, bn in enumerate(vis_params['bname']):
        ant1, ant2 = bn.split('-')
        idx1 = vis_params['antenna_order'].index(int(ant1))
        idx2 = vis_params['antenna_order'].index(int(ant2))
        total_delay[:, bni] = vis_params['baseline_cable_delays'][bni] + ((ant_bw[:, idx1]-ant_bw[:, idx2])*u.m/c.c).to_value(u.nanosecond)

    return total_delay[0]


# input json as argument 1, pointing dec as argument 2 in degrees
f = open(sys.argv[1])
dd = json.load(f)
try:
    tobs = dd[list(dd.keys())[0]]['mjds']
except TypeError:
    tobs = dd['mjds']
pt_dec = float(sys.argv[2])
f.close()

tt = Time(str(tobs),format='mjd')
params = gen_vis_params(tt)
print(tt.isot)
buvw, ant_bw = calculate_uvw_and_geodelay(params, pt_dec*u.deg)
d = get_total_delay(ant_bw,params)
f = open("delays.dat","w")
for dv in d:
    f.write('%f\n'%(dv))
f.close()

               
