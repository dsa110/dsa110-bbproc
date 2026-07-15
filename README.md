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
  Run with and without to compare. v1 implements the SK detector
  (dominant in practice); bandpass/group/sumthreshold are not ported.
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
