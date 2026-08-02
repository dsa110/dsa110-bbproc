#!/usr/bin/env python3
"""IQUV dynamic-spectrum plots + .npy dumps from bbproc toolkit --iquv fils.

Reads <base>_{I,Q,U,V}.fil (float32 SIGPROC, produced by `toolkit --iquv`),
dedisperses to --dm, downsamples in time to ~--tres, bins in frequency by
--fbin, writes the four dynamic spectra (+ freq/time axes) as .npy, and a
2x2 IQUV figure.

Example:
  python3 tools/plot_iquv_dynspec.py \
      --dir /dataz/dsa110/candidates/260715twmx/filterbank \
      --base 260715twmx --dm 1702.610474 --tres 2e-3 --fbin 32
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_fil import read_fil  # noqa: E402


def dedisperse(data, hdr, dm):
    """Integer-sample intra-channel dedispersion, ref = band top (fch1)."""
    nt, nch = data.shape
    f = hdr["fch1"] + np.arange(nch) * hdr["foff"]        # MHz, descending
    fg = f * 1e-3
    dt_ms = 4.15 * dm * (fg ** -2 - (f[0] * 1e-3) ** -2)  # per channel
    sh = np.round(dt_ms * 1e-3 / hdr["tsamp"]).astype(int)
    out = np.empty_like(data)
    for c in range(nch):
        out[:, c] = np.roll(data[:, c], -sh[c])
    return out


def block_mean(a, axis, n):
    """Average non-overlapping blocks of length n along axis (drops remainder)."""
    if n <= 1:
        return a
    a = np.moveaxis(a, axis, 0)
    k = (a.shape[0] // n) * n
    a = a[:k].reshape(k // n, n, *a.shape[1:]).mean(axis=1)
    return np.moveaxis(a, 0, axis)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="directory holding the fils")
    ap.add_argument("--base", required=True, help="event base name (no _I.fil)")
    ap.add_argument("--dm", type=float, required=True, help="DM (pc/cc)")
    ap.add_argument("--tres", type=float, default=2e-3,
                    help="target time resolution [s] (default 2e-3)")
    ap.add_argument("--fbin", type=int, default=32,
                    help="frequency bin factor (channels averaged) [32]")
    ap.add_argument("--twin", type=float, default=0.06,
                    help="plot time window (full width) around burst [s]")
    ap.add_argument("--out", default=None, help="output PNG path")
    args = ap.parse_args()

    stokes = ["I", "Q", "U", "V"]
    cubes = {}
    hdr = None
    for s in stokes:
        p = os.path.join(args.dir, f"{args.base}_{s}.fil")
        h, d = read_fil(p)
        hdr = h
        cubes[s] = dedisperse(d, h, args.dm)
        print(f"loaded {p}: {d.shape}")

    tsamp = hdr["tsamp"]
    tdown = max(1, int(round(args.tres / tsamp)))
    tres = tsamp * tdown
    print(f"time downsample x{tdown} -> {tres*1e3:.3f} ms; "
          f"freq bin x{args.fbin} -> {hdr['nchans']//args.fbin} chans "
          f"({abs(hdr['foff'])*args.fbin:.3f} MHz)")

    ds = {}
    for s in stokes:
        a = block_mean(cubes[s], 0, tdown)      # time
        a = block_mean(a, 1, args.fbin)         # freq
        ds[s] = a

    nt, nf = ds["I"].shape
    freq = hdr["fch1"] + (np.arange(nf) + 0.5) * hdr["foff"] * args.fbin  # MHz
    time = (np.arange(nt) + 0.5) * tres                                   # s

    # per-channel off-burst baseline removal: for each Stokes subtract the
    # median over off-burst time bins per channel. This strips the system
    # bandpass from I and any per-channel offset from Q/U/V, leaving the
    # burst's excess power (the quantity used for polarization analysis).
    prof0 = ds["I"].mean(axis=1)
    pk0 = int(np.argmax(prof0 - np.median(prof0)))
    guard = max(4, int(round(0.02 / tres)))       # ~20 ms guard each side
    offmask = np.ones(nt, bool)
    offmask[max(0, pk0 - guard):pk0 + guard + 1] = False
    for s in stokes:
        ds[s] = ds[s] - np.median(ds[s][offmask], axis=0, keepdims=True)

    # --- save .npy ---
    for s in stokes:
        np.save(os.path.join(args.dir, f"{args.base}_dynspec_{s}.npy"), ds[s])
    np.save(os.path.join(args.dir, f"{args.base}_dynspec_freq_MHz.npy"), freq)
    np.save(os.path.join(args.dir, f"{args.base}_dynspec_time_s.npy"), time)
    print("saved .npy: "
          f"{args.base}_dynspec_{{I,Q,U,V}}.npy (shape {nt}x{nf}), "
          "freq_MHz, time_s")

    # --- locate burst (band-averaged I, off-burst detrended) ---
    prof = ds["I"].mean(axis=1)
    base = np.median(prof)
    pk = int(np.argmax(prof - base))
    t0 = time[pk]
    half = int(round(args.twin / 2 / tres))
    a, b = max(0, pk - half), min(nt, pk + half + 1)
    print(f"burst at bin {pk} (t={t0*1e3:.2f} ms); plotting [{a},{b})")

    # --- figure: 2x2 IQUV dynamic spectra over the burst window ---
    ext = [time[a] * 1e3, time[b - 1] * 1e3, freq[-1], freq[0]]  # ms, MHz
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    for ax, s in zip(axes.ravel(), stokes):
        win = ds[s][a:b].T                      # [freq, time]
        if s == "I":
            vmax = np.nanpercentile(win, 99.5)
            vmin = np.nanpercentile(win, 2.0)
            cmap, kw = "viridis", dict(vmin=vmin, vmax=vmax)
        else:
            v = np.nanpercentile(np.abs(win), 99.0)
            cmap, kw = "RdBu_r", dict(vmin=-v, vmax=v)
        im = ax.imshow(win, aspect="auto", origin="upper", extent=ext,
                       cmap=cmap, interpolation="nearest", **kw)
        ax.set_title(f"Stokes {s}")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    for ax in axes[-1]:
        ax.set_xlabel("time [ms]")
    for ax in axes[:, 0]:
        ax.set_ylabel("freq [MHz]")
    fig.suptitle(f"{args.base}  IQUV dynamic spectra  "
                 f"(DM={args.dm:.2f}, {tres*1e3:.2f} ms, "
                 f"{abs(hdr['foff'])*args.fbin:.2f} MHz)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = args.out or os.path.join(args.dir, f"{args.base}_iquv_dynspec.png")
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
