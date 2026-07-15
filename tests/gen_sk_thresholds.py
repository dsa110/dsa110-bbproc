#!/usr/bin/env python3
"""Regenerate the SK threshold table embedded in src/toolkit.cu.

Reproduces dsa110-rt dsart/rfi/sk.py::_mc_sk_thresholds exactly:
|E|^2 ~ Exp(1), SK = ((M+1)/(M-1)) * (M*S2/S1^2 - 1), 1e6 trials,
two-sided FAR quantiles. Chunked so M=4096 fits in RAM.
"""
import numpy as np

FAR = 1e-4
N = 1_000_000
CHUNK = 50_000

rng = np.random.default_rng(1234)
print("static const SkThresh SK_TABLE[] = {")
for M in (64, 256, 1024, 4096):
    sk = np.empty(N)
    for i in range(0, N, CHUNK):
        p = rng.exponential(1.0, size=(CHUNK, M))
        s1 = p.sum(1)
        s2 = (p * p).sum(1)
        sk[i:i + CHUNK] = (M + 1.0) / (M - 1.0) * (M * s2 / (s1 * s1) - 1.0)
    lo = np.quantile(sk, FAR / 2)
    hi = np.quantile(sk, 1 - FAR / 2)
    print(f"    {{{M}, {lo:.6f}f, {hi:.6f}f}},")
print("};")
