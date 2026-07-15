#!/bin/bash
# GPU round-trip test: fake a dispersed pulse at (l, m), beamform back
# toward it, and check the pulse is recovered — and that pointing away
# kills it. Needs a GPU + ~5 GB scratch. Run from the repo root.
set -eo pipefail

SCRATCH=${BBPROC_TEST_DIR:-/tmp/bbproc_test}
GPU=${BBPROC_TEST_GPU:-0}
EVT=fake0000test
L=0.002
M=-0.001
DM=300
NBLOCKS=4          # 4 blocks = 0.54 s per fragment; 4 subbands below
SBRANGE=0-3        # subband subset keeps the test < 5 GB / < 1 min

mkdir -p "$SCRATCH"
rm -f "$SCRATCH"/${EVT}_sb*_data.out "$SCRATCH"/${EVT}_sb*.json

echo "== faking event (subbands $SBRANGE, $NBLOCKS blocks) =="
./fake_voltages -O "$SCRATCH" -E $EVT --nblocks $NBLOCKS --sb $SBRANGE \
    --l $L --m $M --dm $DM --width 16 --amp 1.2 --t0 0.15 --gpu $GPU

# Synthetic antpos: the faker used antpos=0 (no cal blob), so the (l,m)
# phasor was unity — beamforming ON and OFF position would be identical.
# For a real geometric test we need a cal blob. Build a synthetic one:
# unit gains, antennas on a 60 m-spaced E-W line + N scatter.
python3 - "$SCRATCH/synth_cal.dat" <<'EOF'
import struct, sys
NANT, NC = 96, 48
with open(sys.argv[1], "wb") as f:
    for a in range(NANT):                     # antpos_e
        f.write(struct.pack("<f", (a - NANT/2) * 60.0))
    for a in range(NANT):                     # antpos_n
        f.write(struct.pack("<f", ((a * 37) % 29 - 14) * 25.0))
    for a in range(NANT):                     # unit gains
        for c in range(NC):
            for p in range(2):
                f.write(struct.pack("<ff", 1.0, 0.0))
EOF

echo "== re-faking with geometry =="
./fake_voltages -O "$SCRATCH" -E $EVT --nblocks $NBLOCKS --sb $SBRANGE \
    --l $L --m $M --dm $DM --width 16 --amp 1.2 --t0 0.15 \
    -w "$SCRATCH/synth_cal.dat" --gpu $GPU

echo "== beamform ON position =="
./toolkit -D "$SCRATCH" -E $EVT -P "$SCRATCH/on.fil" --l $L --m $M \
    -w "$SCRATCH/synth_cal.dat" --core all --tscrunch 4 --gpu $GPU

echo "== beamform OFF position =="
./toolkit -D "$SCRATCH" -E $EVT -P "$SCRATCH/off.fil" --l 0 --m 0 \
    -w "$SCRATCH/synth_cal.dat" --core all --tscrunch 4 --gpu $GPU

echo "== beamform ON + dedispersion + RFI flagging =="
./toolkit -D "$SCRATCH" -E $EVT -P "$SCRATCH/on_dm.fil" --l $L --m $M \
    -w "$SCRATCH/synth_cal.dat" --core all --tscrunch 4 --dm $DM --rfi \
    --gpu $GPU

echo "== verifying =="
python3 - "$SCRATCH/on.fil" "$SCRATCH/off.fil" "$SCRATCH/on_dm.fil" <<'EOF'
import struct, sys
import numpy as np

def read_fil(path):
    with open(path, "rb") as f:
        raw = f.read()
    # minimal SIGPROC parse: find HEADER_END, read int keys we need
    def find_int(key):
        i = raw.find(key.encode())
        return struct.unpack("<i", raw[i+len(key):i+len(key)+4])[0]
    end = raw.find(b"HEADER_END") + len(b"HEADER_END")
    nchans = find_int("nchans")
    data = np.frombuffer(raw[end:], dtype=np.float32)
    return data.reshape(-1, nchans)

on, off, ondm = (read_fil(p) for p in sys.argv[1:4])

def pulse_snr(a):
    ts = a.sum(axis=1)
    med = np.median(ts)
    mad = np.median(np.abs(ts - med)) * 1.4826 + 1e-12
    return (ts - med).max() / mad, int(np.argmax(ts))

snr_on, t_on = pulse_snr(on)
snr_off, _ = pulse_snr(off)
snr_dm, t_dm = pulse_snr(ondm)
print(f"ON  pointing: band-summed peak S/N = {snr_on:8.1f} at sample {t_on}")
print(f"OFF pointing: band-summed peak S/N = {snr_off:8.1f}")
print(f"ON + dedisp : band-summed peak S/N = {snr_dm:8.1f} at sample {t_dm}")

ok = True
if snr_dm < 10:
    print("FAIL: dedispersed on-position pulse not recovered"); ok = False
if snr_dm < 1.5 * snr_on:
    print("FAIL: dedispersion did not sharpen the pulse"); ok = False
if snr_on < 3 * snr_off:
    print("FAIL: off-position suppression too weak "
          "(coherent phasing broken?)"); ok = False
# expected arrival: t0=0.15 s / (4 * 32.768us) after dedispersion
exp = int(0.15 / (4 * 32.768e-6))
if abs(t_dm - exp) > 50:
    print(f"FAIL: pulse at sample {t_dm}, expected ~{exp}"); ok = False
print("ROUNDTRIP", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
EOF
