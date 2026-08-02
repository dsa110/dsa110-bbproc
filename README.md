# dsa110-bbproc

Offline baseband processing of **DSA-110 M8 voltage dumps** — the raw
int4 voltage fragments the dsa110-rt `voltage_retention` service stages
on each corr node and C3 collects to h23 under
`/dataz/dsa110/candidates/<event>/Level2/voltages/`.

> The pre-M8 code (legacy T3-era beamformers, splicers, correlators)
> lives unchanged on the **`legacy`** branch.

## Tools

### `toolkit`

One binary, two modes (`./toolkit -h` for everything):

**Coherent filterbank** — the headline mode. Reads a candidate's 16
subband fragments, coherently beamforms the **core antennas** (82 of
96; outriggers excluded, same station≤102 definition as dsa110-rt)
toward a candidate `(l, m)` offset from the meridian phase centre, and
writes one full-band 6144-channel SIGPROC filterbank:

```sh
./toolkit -D /dataz/dsa110/candidates/EVT/Level2/voltages -E EVT \
          -P EVT.fil --l -0.0034 --m 0.0012 \
          -w beamformer_weights_XXX.dat --phase-only \
          --tscrunch 8 [--dm 168.8] [--rfi]
```

* Voltages are **raw off the SNAPs** — this tool applies calibration
  (`-w`, the same 74,496-byte `beamformer_weights_*.dat` blob the
  realtime system distributes), antenna exclusion (`--core`, `-f`,
  zero-gain cal entries), and optionally RFI flagging.
* `(l, m)` phasor convention matches the dsa110-rt injection system
  (`dsart/inject/online.py`): a source at `(l, m)` carries
  `e^{+2πi ν (E·l + N·m)/c}`; we beamform with the conjugate.
* `--rfi` enables Spectral-Kurtosis flagging with the **same statistic
  and thresholds as the realtime flagger** (`dsart/rfi/sk.py`, Nita &
  Gary 2010, MC thresholds at FAR 1e-4, per (ant, ch, pol, window)).
  Run with and without to compare. Only the SK detector is ported
  (dominant in practice); bandpass/group/sumthreshold are not.

  **v2 (2026-08-02) — the v1 flagger made things worse; do not use
  those products.** Two bugs, both now fixed:

  1. v1 rolled the per-channel SK trips up into a per-`(win, ant, pol)`
     verdict and killed the antenna's whole 384-channel subband for the
     8.4 ms window once >1% of its channels tripped. Narrowband RFI in
     4 channels therefore threw away 380 good ones. v2 masks per
     `(win, ant, ch, pol)`, matching the realtime granularity.
  2. v1 renormalized the beam **voltage** by `n_ref/n_live`, which holds
     the coherent signal amplitude fixed and lets the noise power ride
     up as `n_ref²/n_live` — so every flagged window got a *higher*
     noise floor than its neighbours, and the windows with the most RFI
     got boosted the most. It also compared two different antenna
     counts (`n_live` from the core cut, `n_ref` from the weights), so
     a window with **nothing flagged** was still scaled by 0.92 in
     power. v2 counts both from the weights per `(ch, pol)` and scales
     the **power** by `n_ref/n_live` (`--rfi-norm noise`, default), so
     the noise floor is stationary and `--rfi` is a bit-exact no-op
     wherever nothing trips. `--rfi-norm amp` restores the v1 exponent
     for A/B tests; `--rfi-min-live` (0.25) blanks a cell outright
     rather than amplifying what little survives.

  Measured against the v1 products on three 2026-08-02 events
  (`260802unoj`, `260802totk`, `260802gunl`), relative to no flagging:

  | | v1 gain | v1 per-window noise | v1 dedisp S/N | v2 gain | v2 noise | v2 S/N |
  |---|---|---|---|---|---|---|
  | unoj | ×0.916 | ×1.28 | ×0.95 | ×1.000 | ×1.02 | ×0.99 |
  | totk | ×0.927 | ×1.37 | ×1.30 † | ×1.000 | ×1.00 | ×0.99 |
  | gunl | ×0.916 | ×1.33 | ×0.91 | ×0.999 | ×0.99 | ×1.00 |

  † not an improvement — v1 *manufactured* a 5.3σ peak in noise that
  reads 4.1σ unflagged. gunl carries a real 23.1σ burst that v1 pushed
  down to 21.1σ and v2 preserves at 23.0σ.

  The mask and the live-antenna counts are now computed in CUDA
  (`k_skmask` / `k_rfiscale`) instead of copying the 9.4 MB autos array
  back per block and running a 1.2 M-iteration host loop.

  Cross-hand Stokes (`--stokes 1..3`, `--iquv`) additionally lock the
  two pols to the same mask — otherwise the two beams are formed from
  different antenna sets and Q/U/V decorrelate in flagged cells.
* Memory: streamed block-by-block — a full 103 GiB event needs
  ~2.5 GB host RAM and ~600 MB GPU.
* Missing fragments are zero-filled (flagged on stdout).

**Visibilities** (legacy `toolkit_dev` parity, single fragment):
correlate all 4656 baselines with optional time integration `-t`,
per-baseline delay removal `-d`, 8× frequency averaging `-a`, and the
baseline-sum filterbank `-p`/`-g`/`-v`. Output formats unchanged from
the legacy toolkit.

### `fake_voltages`

Synthesizes an M8-format event (fragments + manifests) containing
thermal int4 noise plus a dispersed pulse at a chosen `(l, m)`, DM,
width, and amplitude — the C++/CUDA analogue of the dsa110-rt injection
capability, for end-to-end testing without a telescope:

```sh
./fake_voltages -O /tmp/fake -E fake0000test --nblocks 4 \
    --l 0.002 --m -0.001 --dm 300 --width 16 --amp 1.2 -w cal.dat
./toolkit -D /tmp/fake -E fake0000test -P test.fil --l 0.002 --m -0.001 \
    -w cal.dat --dm 300
```

## Build

On h23 (CUDA 11.1, RTX 2080 Ti):

```sh
make            # toolkit + fake_voltages
make test       # GPU round-trip: fake event -> filterbank -> S/N checks
```

## Format ground truth

Pinned to dsa110-rt (see `src/bbproc.h` for chapter and verse):

* fragment = N × 288 MiB blocks, `[2048 pkt, 96 ant, 384 ch, 2 t, 2 pol]`
  int4-complex (real = low nibble), 32.768 µs native samples;
* band: channel 0 of sb00 = **1498.75 MHz** (system channel 1024),
  descending by 250/8192 MHz across 6144 channels;
* cal blob: `antpos_e[96], antpos_n[96], gains[96][48][2][2]` float32,
  pol order [B, A], each coarse channel spans 8 fine channels.
