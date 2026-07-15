#!/usr/bin/env python3
"""Candidate inspection plot for dsa110-bbproc filterbanks.

Reads a SIGPROC filterbank (as written by ``toolkit -P``) and renders a
single multi-panel PNG highlighting a candidate:

  1. dedispersed waterfall (full file, decimated) with the candidate
     time marked;
  2. dedispersed band-summed time series in robust S/N units;
  3. DM-time "bowtie": S/N vs trial DM around the candidate DM;
  4. zoomed dedispersed waterfall around the pulse;
  5. on-pulse minus off-pulse spectrum;
  6. a text card with the numbers (event, DM, S/N, MJD, ...).

Per-channel baselines are removed (median) and channels scaled by MAD;
dead channels (zero variance) are masked. Dedispersion here is
incoherent integer-shift at the file resolution — fine for inspection
(the toolkit can pre-shift coherently within channels via --dm; pass
--already-dedispersed in that case).

Designed to be called by the C3 integration with values from the C2
trigger row; also usable by hand:

  python3 tools/plot_fil.py EVT.fil --dm 168.8 --t0 1.4 --width 4 \\
      --event 260714spjx --out EVT.png

Dependencies: numpy + matplotlib only (dsart_h23 env on h23 has both).
"""

from __future__ import annotations

import argparse
import struct
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

DISP_MS = 4.15  # ms; dt = DISP_MS * DM * (f_GHz^-2 - fref_GHz^-2)


# ---------------------------------------------------------------------------
# SIGPROC IO
# ---------------------------------------------------------------------------

_INT_KEYS = {"telescope_id", "machine_id", "data_type", "nchans", "nbits",
             "nifs", "barycentric", "pulsarcentric", "nbeams", "ibeam"}
_DBL_KEYS = {"fch1", "foff", "tstart", "tsamp", "az_start", "za_start",
             "src_raj", "src_dej", "refdm", "period"}
_STR_KEYS = {"source_name", "rawdatafile"}


def read_fil(path):
    """Minimal SIGPROC reader -> (header dict, data [ntime, nchan] f32)."""
    hdr = {}
    with open(path, "rb") as f:
        def rstr():
            (n,) = struct.unpack("<i", f.read(4))
            if not 0 < n < 128:
                raise ValueError(f"bad header string length {n}")
            return f.read(n).decode("ascii", "replace")

        if rstr() != "HEADER_START":
            raise ValueError("not a SIGPROC filterbank")
        while True:
            key = rstr()
            if key == "HEADER_END":
                break
            if key in _INT_KEYS:
                (hdr[key],) = struct.unpack("<i", f.read(4))
            elif key in _DBL_KEYS:
                (hdr[key],) = struct.unpack("<d", f.read(8))
            elif key in _STR_KEYS:
                hdr[key] = rstr()
            else:
                raise ValueError(f"unhandled header key {key!r}")
        payload = f.read()

    nchan = hdr["nchans"]
    nbits = hdr.get("nbits", 32)
    dtype = {32: np.float32, 8: np.uint8, 16: np.uint16}.get(nbits)
    if dtype is None:
        raise ValueError(f"nbits={nbits} unsupported")
    data = np.frombuffer(payload, dtype=dtype)
    ntime = data.size // nchan
    return hdr, data[: ntime * nchan].reshape(ntime, nchan).astype(np.float32)


# ---------------------------------------------------------------------------
# processing
# ---------------------------------------------------------------------------

def clean_normalize(data):
    """Per-channel median subtraction + MAD scaling; mask dead channels.

    Returns (norm [ntime, nchan], good_channel_mask)."""
    med = np.median(data, axis=0)
    mad = np.median(np.abs(data - med), axis=0) * 1.4826
    good = mad > 0
    out = np.zeros_like(data)
    out[:, good] = (data[:, good] - med[good]) / mad[good]
    return out, good


def dedisperse(data, freqs_mhz, dm, tsamp_s, fref_mhz):
    """Integer-shift incoherent dedispersion (aligns to fref)."""
    if dm == 0:
        return data
    f_ghz = freqs_mhz * 1e-3
    dt = DISP_MS * 1e-3 * dm * (f_ghz ** -2.0 - (fref_mhz * 1e-3) ** -2.0)
    shifts = np.round(dt / tsamp_s).astype(int)
    out = np.zeros_like(data)
    n = data.shape[0]
    for ch, s in enumerate(shifts):
        if s == 0:
            out[:, ch] = data[:, ch]
        elif 0 < s < n:
            out[: n - s, ch] = data[s:, ch]
        elif -n < s < 0:
            # negative residual shift (bowtie trials below the candidate DM)
            out[-s:, ch] = data[: n + s, ch]
    return out


def snr_series(dd, good, width):
    """Band-summed, boxcar-matched time series in robust S/N units."""
    ts = dd[:, good].sum(axis=1)
    if width > 1:
        kern = np.ones(width, dtype=np.float32)
        ts = np.convolve(ts, kern, mode="same") / np.sqrt(width)
    med = np.median(ts)
    mad = np.median(np.abs(ts - med)) * 1.4826 + 1e-12
    return (ts - med) / mad


def decimate2d(a, max_t=1200, max_f=512):
    """Block-average down to a plottable size."""
    t_fac = max(1, a.shape[0] // max_t)
    f_fac = max(1, a.shape[1] // max_f)
    nt = (a.shape[0] // t_fac) * t_fac
    nf = (a.shape[1] // f_fac) * f_fac
    return a[:nt, :nf].reshape(nt // t_fac, t_fac, nf // f_fac,
                               f_fac).mean(axis=(1, 3)), t_fac, f_fac


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("fil")
    ap.add_argument("--dm", type=float, default=0.0,
                    help="candidate DM (pc/cc)")
    ap.add_argument("--t0", type=float, default=None,
                    help="candidate time, s from file start (dedispersed, "
                         "referenced to the highest channel); default: "
                         "brightest sample")
    ap.add_argument("--width", type=int, default=4,
                    help="candidate boxcar width in .fil samples [4]")
    ap.add_argument("--event", default=None, help="event name for the title")
    ap.add_argument("--lm", nargs=2, type=float, default=None,
                    metavar=("L", "M"), help="beamformed (l,m), for the card")
    ap.add_argument("--out", default=None, help="output PNG [<fil>.png]")
    ap.add_argument("--already-dedispersed", action="store_true",
                    help="the .fil was written with toolkit --dm (skip "
                         "dedispersion here, keep DM for the card/bowtie)")
    ap.add_argument("--dm-trials", type=int, default=64)
    ap.add_argument("--dm-span", type=float, default=0.25,
                    help="bowtie half-span as a fraction of DM [0.25]")
    ap.add_argument("--zoom-widths", type=float, default=64,
                    help="zoom window half-width in units of --width [64]")
    args = ap.parse_args()

    hdr, raw = read_fil(args.fil)
    nchan = hdr["nchans"]
    tsamp = hdr["tsamp"]
    freqs = hdr["fch1"] + np.arange(nchan) * hdr["foff"]   # descending
    fref = freqs.max()
    ntime = raw.shape[0]
    event = args.event or hdr.get("source_name", "candidate")

    norm, good = clean_normalize(raw)
    n_dead = int((~good).sum())

    dd = norm if args.already_dedispersed else dedisperse(
        norm, freqs, args.dm, tsamp, fref)
    ts = snr_series(dd, good, args.width)

    if args.t0 is not None:
        i0 = int(round(args.t0 / tsamp))
        # allow the marker to snap to the local max within +-width*4
        lo = max(0, i0 - 4 * args.width)
        hi = min(ntime, i0 + 4 * args.width + 1)
        ipk = lo + int(np.argmax(ts[lo:hi]))
    else:
        ipk = int(np.argmax(ts))
    snr_pk = float(ts[ipk])
    t_pk = ipk * tsamp

    # ---- bowtie: S/N vs trial DM (on decimated-in-freq data for speed) ----
    dm_trials = None
    bowtie = None
    if args.dm > 0:
        dm_trials = np.linspace(args.dm * (1 - args.dm_span),
                                args.dm * (1 + args.dm_span),
                                args.dm_trials)
        # residual-DM trick: dedisperse the already-dedispersed data by
        # (trial - dm); cheap on a frequency-decimated copy.
        dec_f = max(1, nchan // 256)
        nf = (nchan // dec_f) * dec_f
        small = dd[:, :nf].reshape(ntime, nf // dec_f, dec_f).mean(axis=2)
        gsmall = good[:nf].reshape(nf // dec_f, dec_f).mean(axis=1) > 0.5
        fsmall = freqs[:nf].reshape(nf // dec_f, dec_f).mean(axis=1)
        bowtie = np.empty((args.dm_trials, ntime), dtype=np.float32)
        for k, dmt in enumerate(dm_trials):
            r = dedisperse(small, fsmall, dmt - args.dm, tsamp, fref)
            bowtie[k] = snr_series(r, gsmall, args.width)

    # ---- figure -------------------------------------------------------------
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[1.0, 2.0, 1.4],
                  width_ratios=[2.2, 1.0, 1.0], hspace=0.28, wspace=0.24)

    t_axis_full = np.arange(ntime) * tsamp

    # (2) time series
    ax_ts = fig.add_subplot(gs[0, 0])
    ax_ts.plot(t_axis_full, ts, lw=0.5, color="#2d3436")
    ax_ts.axvline(t_pk, color="#d63031", lw=1.0, alpha=0.8)
    ax_ts.set_xlim(0, t_axis_full[-1])
    ax_ts.set_ylabel("S/N")
    ax_ts.set_title(f"dedispersed time series (boxcar {args.width})")
    ax_ts.annotate(f"S/N {snr_pk:.1f}", (t_pk, snr_pk),
                   textcoords="offset points", xytext=(8, -4),
                   color="#d63031", fontsize=10, fontweight="bold")

    # (1) full waterfall
    ax_wf = fig.add_subplot(gs[1, 0], sharex=ax_ts)
    wf, t_fac, f_fac = decimate2d(dd)
    vmax = np.percentile(wf, 99.5)
    ax_wf.imshow(wf.T, aspect="auto", origin="upper", cmap="viridis",
                 vmin=-1.0, vmax=max(vmax, 1.0),
                 extent=[0, ntime * tsamp, freqs.min(), freqs.max()])
    ax_wf.axvline(t_pk, color="#d63031", lw=0.8, alpha=0.7)
    ax_wf.set_ylabel("frequency (MHz)")
    ax_wf.set_xlabel("time (s)")
    ax_wf.set_title(f"dedispersed waterfall (DM {args.dm:.2f})")

    # (4) zoom waterfall
    ax_zm = fig.add_subplot(gs[1, 1])
    half = int(args.zoom_widths * args.width)
    zlo, zhi = max(0, ipk - half), min(ntime, ipk + half)
    zoom = dd[zlo:zhi]
    zwf, _, _ = decimate2d(zoom, max_t=400, max_f=512)
    ax_zm.imshow(zwf.T, aspect="auto", origin="upper", cmap="viridis",
                 vmin=-1.0, vmax=max(np.percentile(zwf, 99.5), 1.0),
                 extent=[zlo * tsamp, zhi * tsamp, freqs.min(), freqs.max()])
    ax_zm.axvline(t_pk, color="#d63031", lw=0.8, alpha=0.7)
    ax_zm.set_title("zoom (dedispersed)")
    ax_zm.set_xlabel("time (s)")

    # (5) on - off spectrum
    ax_sp = fig.add_subplot(gs[1, 2])
    on = dd[max(0, ipk - args.width): ipk + args.width + 1].mean(axis=0)
    off_sl = np.r_[max(0, ipk - 40 * args.width): max(0, ipk - 10 * args.width)]
    off = dd[off_sl].mean(axis=0) if off_sl.size else np.zeros(nchan)
    spec = np.where(good, on - off, np.nan)
    dec_f = max(1, nchan // 512)
    nf = (nchan // dec_f) * dec_f
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN blocks
        spec_d = np.nanmean(spec[:nf].reshape(nf // dec_f, dec_f), axis=1)
    f_d = freqs[:nf].reshape(nf // dec_f, dec_f).mean(axis=1)
    ax_sp.plot(spec_d, f_d, lw=0.6, color="#0984e3")
    ax_sp.axvline(0, color="#b2bec3", lw=0.5)
    ax_sp.set_ylim(freqs.min(), freqs.max())
    ax_sp.set_title("on-off pulse spectrum")
    ax_sp.set_xlabel("S/N per channel")

    # (3) bowtie
    ax_bt = fig.add_subplot(gs[2, 0], sharex=ax_ts)
    if bowtie is not None:
        ax_bt.imshow(bowtie, aspect="auto", origin="lower", cmap="magma",
                     vmin=0, vmax=max(np.percentile(bowtie, 99.9), 5.0),
                     extent=[0, ntime * tsamp, dm_trials[0], dm_trials[-1]])
        ax_bt.axhline(args.dm, color="#00cec9", lw=0.7, alpha=0.8)
        ax_bt.axvline(t_pk, color="#d63031", lw=0.8, alpha=0.7)
        ax_bt.set_ylabel("trial DM (pc cm$^{-3}$)")
        ax_bt.set_title("DM-time")
    else:
        ax_bt.text(0.5, 0.5, "no DM given — bowtie skipped",
                   transform=ax_bt.transAxes, ha="center", color="#636e72")
    ax_bt.set_xlabel("time (s)")

    # (3b) bowtie profile: S/N at t_pk vs DM
    ax_bp = fig.add_subplot(gs[2, 1])
    if bowtie is not None:
        ax_bp.plot(dm_trials, bowtie[:, ipk], color="#6c5ce7", lw=1.0)
        ax_bp.axvline(args.dm, color="#00cec9", lw=0.8)
        ax_bp.set_xlabel("trial DM (pc cm$^{-3}$)")
        ax_bp.set_ylabel("S/N at candidate time")
        ax_bp.set_title("DM profile")
    else:
        ax_bp.axis("off")

    # (6) info card
    ax_tx = fig.add_subplot(gs[0, 1:])
    ax_tx.axis("off")
    lines = [
        f"event        : {event}",
        f"file         : {args.fil.split('/')[-1]}",
        f"DM           : {args.dm:.3f} pc cm$^{{-3}}$"
        + ("  (pre-dedispersed .fil)" if args.already_dedispersed else ""),
        f"peak S/N     : {snr_pk:.1f}  (boxcar {args.width} x "
        f"{tsamp * 1e6:.0f} us)",
        f"t_peak       : {t_pk:.4f} s from file start",
        f"MJD start    : {hdr.get('tstart', 0):.9f}",
        f"band         : {freqs.max():.2f} - {freqs.min():.2f} MHz, "
        f"{nchan} ch",
        f"dead channels: {n_dead}",
    ]
    if args.lm:
        lines.insert(2, f"(l, m)       : ({args.lm[0]:+.5f}, "
                        f"{args.lm[1]:+.5f}) rad")
    ax_tx.text(0.02, 0.95, "\n".join(lines), transform=ax_tx.transAxes,
               va="top", family="monospace", fontsize=11)

    # bowtie bottom-right: leave for future (polarisation etc.)
    ax_sp2 = fig.add_subplot(gs[2, 2])
    prof = ts[zlo:zhi]
    ax_sp2.plot(np.arange(zlo, zhi) * tsamp, prof, lw=0.7, color="#2d3436")
    ax_sp2.axvline(t_pk, color="#d63031", lw=0.8, alpha=0.8)
    ax_sp2.set_title("pulse profile (zoom)")
    ax_sp2.set_xlabel("time (s)")
    ax_sp2.set_ylabel("S/N")

    fig.suptitle(f"{event} — dsa110-bbproc candidate inspection",
                 fontsize=14, fontweight="bold")

    out = args.out or (args.fil.rsplit(".", 1)[0] + ".png")
    fig.savefig(out, dpi=110, bbox_inches="tight")
    print(f"wrote {out}  (peak S/N {snr_pk:.1f} at {t_pk:.4f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
