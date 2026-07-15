// -*- c++ -*-
// dsa110-bbproc toolkit — offline processing of DSA-110 M8 voltage dumps.
//
// Two modes, one binary:
//
// 1. COHERENT FILTERBANK (-P out.fil): reads a candidate event's voltage
//    fragments (all 16 subbands, or any subset present), coherently
//    beamforms the CORE antennas toward a candidate (l, m) offset from
//    the meridian phase centre (same phasor convention as the dsa110-rt
//    injection system), applies the beamformer-weights cal blob, and
//    writes one full-band (6144-channel) SIGPROC filterbank.
//    Optional: realtime-equivalent SK RFI flagging (--rfi), intra-channel
//    integer-shift dedispersion (--dm), time integration (--tscrunch).
//
// 2. VISIBILITIES (-o out.vis, legacy toolkit_dev parity): correlates a
//    single fragment (all 4656 baselines), with optional time
//    integration (-t), per-baseline delays (-d), 8x frequency averaging
//    (-a), and the legacy baseline-sum filterbank (-p / -g / -v).
//
// Optimization vs the legacy toolkit_dev.cu: beamforming is a direct
// O(NANT) coherent sum (not an O(NBASE)=4656-baseline correlation), the
// input is streamed block-by-block (288 MiB device footprint instead of
// a 4.5 GB monolithic buffer), and unpack+weight+sum are fused in one
// kernel. A full 16-fragment (103 GiB) event fits any GPU.
//
// Voltages are RAW off the SNAPs: uncalibrated, unflagged, no injections.
// Everything downstream of fada (cal, flags, RFI) is this tool's job.

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <getopt.h>
#include <unistd.h>

#include <string>
#include <vector>

#include "bbproc.h"

using namespace bbproc;

// ---------------------------------------------------------------------------
// SK RFI flagging (mirrors dsa110-rt dsart/rfi/sk.py)
//
//   SK_M = ((M+1)/(M-1)) * (M * S2 / S1^2 - 1)   per (ant, ch, pol, window)
//
// Thresholds are the Nita & Gary (2010) confidence interval at two-sided
// FAR = 1e-4, computed by the SAME Monte-Carlo model dsart uses
// (|E|^2 ~ Exp(1), 1e6 trials; dsart/rfi/sk.py::_mc_sk_thresholds).
// Regenerate with tests/gen_sk_thresholds.py.
// ---------------------------------------------------------------------------
struct SkThresh { int m; float lo, hi; };
static const SkThresh SK_TABLE[] = {
    {64,   0.415286f, 3.010445f},
    {256,  0.639271f, 1.794989f},
    {1024, 0.792937f, 1.308578f},
    {4096, 0.889956f, 1.133969f},
};

static bool sk_thresholds(int m, float *lo, float *hi) {
    for (const auto &t : SK_TABLE)
        if (t.m == m) { *lo = t.lo; *hi = t.hi; return true; }
    return false;
}

// ---------------------------------------------------------------------------
// Device kernels — filterbank path
// ---------------------------------------------------------------------------

// Per-(ant, ch, pol) auto-power moments over M-sample windows of one block.
// grid: (NANT, nwin)   block: NCHAN_SB threads (one per channel)
// out: [nwin][NANT][NCHAN_SB][NPOL][2] (S1, S2) float
__global__ void k_autos(const unsigned char *__restrict__ raw,
                        float *__restrict__ out, int m_win) {
    const int ant = blockIdx.x;
    const int win = blockIdx.y;
    const int ch = threadIdx.x;
    float s1[NPOL] = {0.f, 0.f}, s2[NPOL] = {0.f, 0.f};
    const int t0 = win * m_win;
    for (int t = t0; t < t0 + m_win; t++) {
        const int pkt = t >> 1, tt = t & 1;
        const size_t base =
            ((((size_t)pkt * NANT + ant) * NCHAN_SB + ch) * NT_PER_PKT + tt) *
            NPOL;
        for (int p = 0; p < NPOL; p++) {
            const unsigned char b = raw[base + p];
            const float re = (float)((char)((b & 0x0F) << 4) >> 4);
            const float im = (float)((char)(b & 0xF0) >> 4);
            const float pw = re * re + im * im;
            s1[p] += pw;
            s2[p] += pw * pw;
        }
    }
    const size_t o =
        ((((size_t)win * NANT + ant) * NCHAN_SB + ch) * NPOL) * 2;
    for (int p = 0; p < NPOL; p++) {
        out[o + 2 * p] = s1[p];
        out[o + 2 * p + 1] = s2[p];
    }
}

// Coherent beamform of one block.
// grid: T_PER_BLOCK   block: NCHAN_SB threads (one per channel)
// w:    [NANT][NPOL][NCHAN_SB][2]  conj(cal x phasor) weights, 0 = excluded
// mask: [nwin][NANT][NPOL] uint8 (1 = flagged) or nullptr
// norm: [nwin] not needed — normalization folded into host scale per window?
//       Normalization: per (win) the number of live antennas varies per
//       (ant,pol) mask; we renormalize by live-antenna count per (win,pol)
//       computed on host into wnorm[nwin][NPOL].
// out:  [T_PER_BLOCK][NCHAN_SB] float, the requested Stokes parameter.
__global__ void k_beamform(const unsigned char *__restrict__ raw,
                           const float *__restrict__ w,
                           const unsigned char *__restrict__ mask,
                           const float *__restrict__ wnorm, int m_win,
                           int stokes, float *__restrict__ out) {
    const int t = blockIdx.x;
    const int ch = threadIdx.x;
    const int pkt = t >> 1, tt = t & 1;
    const int win = t / m_win;

    float accr[NPOL] = {0.f, 0.f}, acci[NPOL] = {0.f, 0.f};
    for (int a = 0; a < NANT; a++) {
        const size_t base =
            ((((size_t)pkt * NANT + a) * NCHAN_SB + ch) * NT_PER_PKT + tt) *
            NPOL;
        for (int p = 0; p < NPOL; p++) {
            const float wr = w[(((size_t)a * NPOL + p) * NCHAN_SB + ch) * 2];
            const float wi =
                w[(((size_t)a * NPOL + p) * NCHAN_SB + ch) * 2 + 1];
            if (wr == 0.f && wi == 0.f) continue;
            if (mask && mask[((size_t)win * NANT + a) * NPOL + p]) continue;
            const unsigned char b = raw[base + p];
            const float re = (float)((char)((b & 0x0F) << 4) >> 4);
            const float im = (float)((char)(b & 0xF0) >> 4);
            accr[p] += re * wr - im * wi;
            acci[p] += re * wi + im * wr;
        }
    }
    const float n0 = wnorm ? wnorm[(size_t)win * NPOL] : 1.f;
    const float n1 = wnorm ? wnorm[(size_t)win * NPOL + 1] : 1.f;
    const float b_r = accr[0] * n0, b_i = acci[0] * n0;
    const float a_r = accr[1] * n1, a_i = acci[1] * n1;
    float v;
    switch (stokes) {
        case 1:  v = (b_r * b_r + b_i * b_i) - (a_r * a_r + a_i * a_i); break;
        case 2:  v = 2.f * (b_r * a_r + b_i * a_i); break;
        case 3:  v = 2.f * (b_i * a_r - b_r * a_i); break;
        default: v = (b_r * b_r + b_i * b_i) + (a_r * a_r + a_i * a_i); break;
    }
    out[(size_t)t * NCHAN_SB + ch] = v;
}

// ---------------------------------------------------------------------------
// Device kernels — visibilities path (legacy toolkit_dev parity)
// ---------------------------------------------------------------------------
constexpr int NBASE = NANT * (NANT + 1) / 2;  // 4656
constexpr int NPTR = 8;                       // 2t x 2pol x r/i floats
constexpr int AV = 8;                         // -a frequency averaging

// promote one packet's int4 bytes to floats [ant][ch][2t x 2pol x ri]
__global__ void k_promote(const unsigned char *__restrict__ input,
                          float *__restrict__ output) {
    const size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned char b = input[i];
    output[2 * i] = (float)((char)((b & 0x0F) << 4) >> 4);
    output[2 * i + 1] = (float)((char)(b & 0xF0) >> 4);
}

// legacy correlator: one packet (2 time samples), all baselines/chans.
// run with NBASE*NCHAN_SB/32 blocks of 32 threads.
__global__ void k_correlate(const float *__restrict__ input,
                            float *__restrict__ output,
                            const int *__restrict__ a1,
                            const int *__restrict__ a2, float scfac,
                            const float *__restrict__ weights) {
    const int iidx = blockIdx.x * 32 + threadIdx.x;
    const int basel = iidx / NCHAN_SB;
    const int ch = iidx % NCHAN_SB;
    if (basel >= NBASE) return;

    const int i1 = a1[basel] * NCHAN_SB * NPTR + ch * NPTR;
    const int i2 = a2[basel] * NCHAN_SB * NPTR + ch * NPTR;

    // cal weights are per 8-fine-channel coarse block: [ant][48][pol][ri]
    float w1[4], w2[4];
    for (int i = 0; i < 4; i++) {
        w1[i] = weights[a1[basel] * (NCAL_COARSE * 4) + (ch / FINE_PER_COARSE) * 4 + i];
        w2[i] = weights[a2[basel] * (NCAL_COARSE * 4) + (ch / FINE_PER_COARSE) * 4 + i];
    }

    float out_t[2][8];
    for (int ti = 0; ti < 2; ti++) {
        int ii = 0;
        for (int p1 = 0; p1 < 2; p1++) {
            for (int p2 = 0; p2 < 2; p2++) {
                const float a1r = input[i1 + ti * 4 + p1 * 2];
                const float a1i = input[i1 + ti * 4 + p1 * 2 + 1];
                const float a2r = input[i2 + ti * 4 + p2 * 2];
                const float a2i = input[i2 + ti * 4 + p2 * 2 + 1];
                const float w1r = a1r * w1[2 * p1] - a1i * w1[2 * p1 + 1];
                const float w1i = a1r * w1[2 * p1 + 1] + a1i * w1[2 * p1];
                const float w2r = a2r * w2[2 * p2] - a2i * w2[2 * p2 + 1];
                const float w2i = a2r * w2[2 * p2 + 1] + a2i * w2[2 * p2];
                out_t[ti][2 * ii] = w1r * w2r + w1i * w2i;
                out_t[ti][2 * ii + 1] = w1r * w2i - w1i * w2r;
                ii++;
            }
        }
    }
    size_t o = ((size_t)basel * NCHAN_SB + ch) * 8;
    for (int i = 0; i < 8; i++) output[o + i] += out_t[0][i] * scfac;
    o += (size_t)NBASE * NCHAN_SB * 8;
    for (int i = 0; i < 8; i++) output[o + i] += out_t[1][i] * scfac;
}

__global__ void k_adder(const float *__restrict__ in, float *__restrict__ out) {
    const size_t i = (size_t)blockIdx.x * 32 + threadIdx.x;
    out[i] = in[i] + in[(size_t)NBASE * NCHAN_SB * 8 + i];
}

__global__ void k_zero(float *buf) {
    buf[(size_t)blockIdx.x * 32 + threadIdx.x] = 0.f;
}

// remove per-baseline delays: multiply by exp(-2 pi i nu tau)
__global__ void k_delay(float *__restrict__ vis,
                        const float *__restrict__ freqs,
                        const float *__restrict__ delays) {
    const int iidx = blockIdx.x * 32 + threadIdx.x;
    const int bci = iidx / 4;
    const int basel = bci / NCHAN_SB;
    const int ch = bci % NCHAN_SB;
    const float arg = -2.f * (float)M_PI * freqs[ch] * delays[basel] * 1e-9f;
    const float c = cosf(arg), s = sinf(arg);
    const float vr = vis[2 * iidx] * c - vis[2 * iidx + 1] * s;
    const float vi = vis[2 * iidx] * s + vis[2 * iidx + 1] * c;
    vis[2 * iidx] = vr;
    vis[2 * iidx + 1] = vi;
}

// 8x frequency averaging of XX and YY
__global__ void k_fscrunch(const float *__restrict__ in,
                           float *__restrict__ out) {
    const int iidx = blockIdx.x * 32 + threadIdx.x;
    const int bcli = iidx / 4;
    const int poli = iidx % 4;
    const int basel = bcli / (NCHAN_SB / AV);
    const int lch = bcli % (NCHAN_SB / AV);
    const int sumss[4] = {0, 1, 6, 7};
    float acc = 0.f;
    for (int i = 0; i < AV; i++)
        acc += in[((size_t)basel * NCHAN_SB + (AV * lch + i)) * 8 + sumss[poli]];
    out[iidx] = acc;
}

// baseline-sum "philterbank": one Stokes value per channel
__global__ void k_reduce(const float *__restrict__ in, float *__restrict__ out,
                         float scfac, const int *__restrict__ a1,
                         const int *__restrict__ a2, int stokes,
                         const float *__restrict__ antpos_e, float minBase) {
    const int ch = blockIdx.x;
    const int tidx = threadIdx.x;
    volatile __shared__ float summer[1024];
    summer[tidx] = 0.f;
    if (tidx < 582) {
        for (int k = 0; k < 8; k++) {
            const int b = tidx + k * 582;
            if (b >= NBASE) continue;
            if (a1[b] == a2[b]) continue;
            if (fabsf(antpos_e[a2[b]] - antpos_e[a1[b]]) <= minBase) continue;
            const size_t o = ((size_t)b * NCHAN_SB + ch) * 8;
            switch (stokes) {
                case 1: summer[tidx] += in[o] - in[o + 6]; break;
                case 2: summer[tidx] += in[o + 2] + in[o + 4]; break;
                case 3: summer[tidx] += in[o + 3] - in[o + 5]; break;
                default: summer[tidx] += in[o] + in[o + 6]; break;
            }
        }
    }
    __syncthreads();
    for (int s = 512; s > 0; s >>= 1) {
        if (tidx < s) summer[tidx] += summer[tidx + s];
        __syncthreads();
    }
    if (tidx == 0) out[ch] = summer[0] * scfac;
}

// ---------------------------------------------------------------------------
// Host: weights
// ---------------------------------------------------------------------------

struct WeightOpts {
    bool have_cal = false;
    CalBlob cal;
    bool phase_only = false;
    bool swap_pol = false;
    std::vector<bool> use;  // per antenna, post core/flag cuts
    double l = 0.0, m = 0.0;
};

// Fill w[NANT][NPOL][NCHAN_SB][2] = conj(cal gain applied) * conj(phasor).
// Excluded antennas/pols get exact zeros (kernel skip).
static void build_bf_weights(const WeightOpts &o, int sb,
                             std::vector<float> *w_out, int *n_used) {
    std::vector<float> &w = *w_out;
    w.assign((size_t)NANT * NPOL * NCHAN_SB * 2, 0.f);
    const double n_term =
        sqrt(fmax(0.0, 1.0 - o.l * o.l - o.m * o.m)) - 1.0;
    int used = 0;
    for (int a = 0; a < NANT; a++) {
        if (!o.use[a]) continue;
        bool ant_counted = false;
        const double E = o.have_cal ? o.cal.antpos_e[a] : 0.0;
        const double N = o.have_cal ? o.cal.antpos_n[a] : 0.0;
        for (int ch = 0; ch < NCHAN_SB; ch++) {
            const double f_hz = freq_mhz(sb, ch) * 1e6;
            // dsa110-rt injection convention: source at (l,m) appears with
            // e^{+2 pi i nu (E l + N m + U n)/c}; beamform with conjugate.
            const double geo =
                2.0 * M_PI * f_hz * (E * o.l + N * o.m + 0.0 * n_term) / CVAC;
            const double pr = cos(geo), pi = -sin(geo);  // conj(phasor)
            for (int p = 0; p < NPOL; p++) {
                double gr = 1.0, gi = 0.0;
                if (o.have_cal) {
                    const int cp = o.swap_pol ? (1 - p) : p;
                    gr = o.cal.gains[a][ch / FINE_PER_COARSE][cp][0];
                    gi = o.cal.gains[a][ch / FINE_PER_COARSE][cp][1];
                    const double g2 = gr * gr + gi * gi;
                    if (g2 == 0.0) continue;  // cal-flagged antenna/pol
                    const double den = o.phase_only ? sqrt(g2) : g2;
                    // conj(g)/den
                    const double tr = gr / den, ti = -gi / den;
                    gr = tr;
                    gi = ti;
                }
                const size_t idx = (((size_t)a * NPOL + p) * NCHAN_SB + ch) * 2;
                w[idx] = (float)(gr * pr - gi * pi);
                w[idx + 1] = (float)(gr * pi + gi * pr);
                if (!ant_counted && (w[idx] != 0.f || w[idx + 1] != 0.f)) {
                    ant_counted = true;
                    used++;
                }
            }
        }
    }
    *n_used = used;
}

// ---------------------------------------------------------------------------
// Host: filterbank mode
// ---------------------------------------------------------------------------

struct FilOpts {
    std::string event_dir, event, out_fil, cal_path, core_path, flag_path;
    std::string single_frag;  // -i: single-fragment mode
    int single_sb = -1;
    double l = 0.0, m = 0.0;
    int tscrunch = 8;
    double dm = 0.0;
    int stokes = 0;
    bool rfi = false;
    int rfi_m = 256;
    bool phase_only = false, swap_pol = false;
    double mjd_override = -1.0;
    int gpu = 0;
    int telescope_id = 0;
};

static int infer_sb_from_name(const std::string &path) {
    size_t p = path.rfind("_sb");
    if (p == std::string::npos || p + 5 > path.size()) return -1;
    return atoi(path.substr(p + 3, 2).c_str());
}

// The realtime system distributes PER-SUBBAND cal blobs
// (beamformer_weights_sb<NN>_<isot>.dat under .../beamformer_weights/
// applied/). If the -w path contains an "sb<NN>" token, substitute the
// current subband; otherwise use the same blob for every subband.
static std::string cal_path_for_sb(const std::string &path, int sb) {
    size_t p = path.rfind("sb");
    if (p != std::string::npos && p + 4 <= path.size() &&
        isdigit(path[p + 2]) && isdigit(path[p + 3])) {
        std::string out = path;
        char d[3];
        snprintf(d, sizeof d, "%02d", sb);
        out[p + 2] = d[0];
        out[p + 3] = d[1];
        return out;
    }
    return path;
}

static int run_filterbank(const FilOpts &opt) {
    CUDA_CHECK(cudaSetDevice(opt.gpu));

    // ---- antenna selection -------------------------------------------------
    WeightOpts wo;
    wo.l = opt.l;
    wo.m = opt.m;
    wo.phase_only = opt.phase_only;
    wo.swap_pol = opt.swap_pol;
    wo.use.assign(NANT, false);
    if (opt.core_path == "all") {
        for (int a = 0; a < NANT; a++) wo.use[a] = true;
    } else {
        std::vector<int> core;
        if (load_int_list(opt.core_path.c_str(), &core)) return 1;
        for (int a : core)
            if (a >= 0 && a < NANT) wo.use[a] = true;
        printf("core antennas: %zu from %s\n", core.size(),
               opt.core_path.c_str());
    }
    if (!opt.flag_path.empty()) {
        std::vector<int> fl;
        if (load_int_list(opt.flag_path.c_str(), &fl)) return 1;
        for (int a : fl)
            if (a >= 0 && a < NANT) wo.use[a] = false;
        printf("flagged %zu antennas from %s\n", fl.size(),
               opt.flag_path.c_str());
    }
    if (!opt.cal_path.empty()) {
        // Load sb00's (or the generic) blob up front for an early error;
        // per-subband blobs are (re)loaded inside the subband loop.
        if (load_cal_blob(cal_path_for_sb(opt.cal_path, 0).c_str(), &wo.cal))
            return 1;
        wo.have_cal = true;
        printf("cal blob: %s (%s)\n", opt.cal_path.c_str(),
               opt.phase_only ? "phase-only" : "amplitude+phase");
    } else {
        printf("WARNING: no cal blob (-w): unit gains, antpos unknown -> "
               "(l,m) phasing DISABLED (zero baselines)\n");
    }

    // ---- fragment inventory ------------------------------------------------
    std::vector<std::string> frags(NSB);
    std::vector<long long> fsizes(NSB, -1);
    int sb_lo = 0, sb_hi = NSB - 1;
    if (!opt.single_frag.empty()) {
        int sb = opt.single_sb >= 0 ? opt.single_sb
                                    : infer_sb_from_name(opt.single_frag);
        if (sb < 0) {
            fprintf(stderr, "cannot infer subband of %s; pass --sb\n",
                    opt.single_frag.c_str());
            return 1;
        }
        sb_lo = sb_hi = sb;
        frags[sb] = opt.single_frag;
        fsizes[sb] = file_size(opt.single_frag.c_str());
    } else {
        for (int sb = 0; sb < NSB; sb++) {
            frags[sb] = frag_path(opt.event_dir, opt.event, sb);
            fsizes[sb] = file_size(frags[sb].c_str());
        }
    }
    long long max_bytes = 0;
    int n_present = 0;
    for (int sb = sb_lo; sb <= sb_hi; sb++) {
        if (fsizes[sb] > 0) {
            n_present++;
            if (fsizes[sb] % BLOCK_BYTES)
                fprintf(stderr, "WARNING: %s size %lld not a whole number of "
                        "blocks\n", frags[sb].c_str(), fsizes[sb]);
            if (fsizes[sb] > max_bytes) max_bytes = fsizes[sb];
        } else {
            printf("subband %02d: missing (%s) -> zeros\n", sb,
                   frags[sb].c_str());
        }
    }
    if (!n_present) {
        fprintf(stderr, "no fragments found\n");
        return 1;
    }
    const int nblocks = (int)(max_bytes / BLOCK_BYTES);
    const long long ntime = (long long)nblocks * T_PER_BLOCK;
    const int nchan_out =
        (sb_lo == sb_hi) ? NCHAN_SB : NCHAN_FULL;
    printf("%d/%d fragments, %d blocks = %lld samples (%.2f s), %d chans\n",
           n_present, sb_hi - sb_lo + 1, nblocks, ntime,
           ntime * TSAMP_S, nchan_out);

    // ---- tstart from manifest ----------------------------------------------
    double tstart_mjd = opt.mjd_override;
    if (tstart_mjd < 0 && !opt.event_dir.empty()) {
        for (int sb = sb_lo; sb <= sb_hi && tstart_mjd < 0; sb++) {
            std::string txt;
            if (!read_text_file(manifest_path(opt.event_dir, opt.event, sb)
                                    .c_str(), &txt))
                continue;
            double mjd = 0, tgt = 0, first = 0;
            if (json_get_double(txt, "mjd_target", &mjd) && mjd > 40000.0 &&
                json_get_double(txt, "target_block_n", &tgt) &&
                json_get_double(txt, "block_n_first", &first)) {
                tstart_mjd =
                    mjd - (tgt - first) * T_PER_BLOCK * TSAMP_S / 86400.0;
            }
        }
    }
    if (tstart_mjd < 0) {
        printf("WARNING: no usable manifest mjd_target; tstart=0 "
               "(pass --mjd to set)\n");
        tstart_mjd = 0.0;
    }

    // ---- SK thresholds -----------------------------------------------------
    float sk_lo = 0.f, sk_hi = 0.f;
    int nwin = 0;
    if (opt.rfi) {
        if (!sk_thresholds(opt.rfi_m, &sk_lo, &sk_hi)) {
            fprintf(stderr, "no SK threshold table entry for M=%d "
                    "(have 64,256,1024,4096)\n", opt.rfi_m);
            return 1;
        }
        if (T_PER_BLOCK % opt.rfi_m) {
            fprintf(stderr, "--rfi-m must divide %d\n", T_PER_BLOCK);
            return 1;
        }
        nwin = T_PER_BLOCK / opt.rfi_m;
        printf("RFI: SK per (ant,ch,pol) M=%d FAR=1e-4 -> [%.4f, %.4f]\n",
               opt.rfi_m, sk_lo, sk_hi);
    }

    // ---- buffers -----------------------------------------------------------
    unsigned char *h_block, *d_raw;
    float *d_w, *d_pow, *d_autos = nullptr, *d_wnorm = nullptr;
    unsigned char *d_mask = nullptr;
    CUDA_CHECK(cudaHostAlloc(&h_block, BLOCK_BYTES, cudaHostAllocDefault));
    CUDA_CHECK(cudaMalloc(&d_raw, BLOCK_BYTES));
    CUDA_CHECK(cudaMalloc(&d_w, (size_t)NANT * NPOL * NCHAN_SB * 2 * 4));
    CUDA_CHECK(cudaMalloc(&d_pow, (size_t)T_PER_BLOCK * NCHAN_SB * 4));
    std::vector<float> h_pow((size_t)T_PER_BLOCK * NCHAN_SB);
    std::vector<float> h_autos;
    std::vector<unsigned char> h_mask;
    std::vector<float> h_wnorm;
    if (opt.rfi) {
        h_autos.resize((size_t)nwin * NANT * NCHAN_SB * NPOL * 2);
        h_mask.resize((size_t)nwin * NANT * NPOL);
        h_wnorm.resize((size_t)nwin * NPOL);
        CUDA_CHECK(cudaMalloc(&d_autos, h_autos.size() * 4));
        CUDA_CHECK(cudaMalloc(&d_mask, h_mask.size()));
        CUDA_CHECK(cudaMalloc(&d_wnorm, h_wnorm.size() * 4));
    }

    // full-band accumulation buffer [ntime][nchan_out] (float32).
    // 23 blocks x 6144 ch = 2.3 GB — fine on h23 host RAM.
    std::vector<float> full((size_t)ntime * nchan_out, 0.f);

    long long tot_flagged_cells = 0, tot_cells = 0;

    // ---- per-subband streaming ---------------------------------------------
    for (int sb = sb_lo; sb <= sb_hi; sb++) {
        if (fsizes[sb] <= 0) continue;
        const int ch_off = (sb_lo == sb_hi) ? 0 : sb * NCHAN_SB;

        if (wo.have_cal && sb != 0) {
            const std::string cp = cal_path_for_sb(opt.cal_path, sb);
            if (load_cal_blob(cp.c_str(), &wo.cal)) {
                fprintf(stderr, "WARNING: %s missing — using previous "
                        "subband's blob\n", cp.c_str());
            }
        }

        std::vector<float> w;
        int n_used = 0;
        build_bf_weights(wo, sb, &w, &n_used);
        if (sb == sb_lo)
            printf("beamforming with %d antennas toward (l,m)=(%.6g,%.6g) rad\n",
                   n_used, opt.l, opt.m);
        CUDA_CHECK(cudaMemcpy(d_w, w.data(), w.size() * 4,
                              cudaMemcpyHostToDevice));

        FILE *f = fopen(frags[sb].c_str(), "rb");
        if (!f) continue;
        const int nb = (int)(fsizes[sb] / BLOCK_BYTES);
        for (int b = 0; b < nb; b++) {
            if (fread(h_block, 1, BLOCK_BYTES, f) != (size_t)BLOCK_BYTES) {
                fprintf(stderr, "short read %s block %d\n", frags[sb].c_str(), b);
                break;
            }
            CUDA_CHECK(cudaMemcpy(d_raw, h_block, BLOCK_BYTES,
                                  cudaMemcpyHostToDevice));

            if (opt.rfi) {
                dim3 g(NANT, nwin);
                k_autos<<<g, NCHAN_SB>>>(d_raw, d_autos, opt.rfi_m);
                CUDA_CHECK(cudaMemcpy(h_autos.data(), d_autos,
                                      h_autos.size() * 4,
                                      cudaMemcpyDeviceToHost));
                // SK mask on host: flag (win, ant, pol) if ANY channel's SK
                // trips (channel-resolved masking costs a [win][ant][ch][pol]
                // mask; realtime combines detectors per (ant,pol) similarly
                // before zeroing — see dsart/rfi/combine.py). We flag the
                // (win, ant, pol) cell when >1% of its channels trip, which
                // suppresses broadband bursts without killing an antenna for
                // a single hot channel.
                const float m_ = (float)opt.rfi_m;
                const float c1 = (m_ + 1.f) / (m_ - 1.f);
                for (int win = 0; win < nwin; win++) {
                    for (int a = 0; a < NANT; a++) {
                        for (int p = 0; p < NPOL; p++) {
                            int ntrip = 0;
                            for (int ch = 0; ch < NCHAN_SB; ch++) {
                                const size_t o =
                                    ((((size_t)win * NANT + a) * NCHAN_SB + ch) *
                                     NPOL + p) * 2;
                                const float s1 = h_autos[o];
                                const float s2 = h_autos[o + 1];
                                if (s1 <= 0.f) continue;
                                const float sk =
                                    c1 * (m_ * s2 / (s1 * s1) - 1.f);
                                if (sk < sk_lo || sk > sk_hi) ntrip++;
                            }
                            const bool bad = ntrip > NCHAN_SB / 100;
                            h_mask[((size_t)win * NANT + a) * NPOL + p] =
                                bad ? 1 : 0;
                            if (bad) tot_flagged_cells++;
                            tot_cells++;
                        }
                    }
                }
                // per-(win,pol) renormalization by live antenna count
                for (int win = 0; win < nwin; win++) {
                    for (int p = 0; p < NPOL; p++) {
                        int live = 0;
                        for (int a = 0; a < NANT; a++)
                            if (wo.use[a] &&
                                !h_mask[((size_t)win * NANT + a) * NPOL + p])
                                live++;
                        h_wnorm[(size_t)win * NPOL + p] =
                            live > 0 ? (float)n_used / (float)live : 0.f;
                    }
                }
                CUDA_CHECK(cudaMemcpy(d_mask, h_mask.data(), h_mask.size(),
                                      cudaMemcpyHostToDevice));
                CUDA_CHECK(cudaMemcpy(d_wnorm, h_wnorm.data(),
                                      h_wnorm.size() * 4,
                                      cudaMemcpyHostToDevice));
            }

            k_beamform<<<T_PER_BLOCK, NCHAN_SB>>>(
                d_raw, d_w, opt.rfi ? d_mask : nullptr,
                opt.rfi ? d_wnorm : nullptr, opt.rfi ? opt.rfi_m : T_PER_BLOCK,
                opt.stokes, d_pow);
            CUDA_CHECK(cudaMemcpy(h_pow.data(), d_pow, h_pow.size() * 4,
                                  cudaMemcpyDeviceToHost));

            // scatter into the full-band buffer
            const long long t0 = (long long)b * T_PER_BLOCK;
            for (int t = 0; t < T_PER_BLOCK; t++) {
                float *dst = &full[(t0 + t) * nchan_out + ch_off];
                memcpy(dst, &h_pow[(size_t)t * NCHAN_SB], NCHAN_SB * 4);
            }
        }
        fclose(f);
        printf("subband %02d done (%s)\n", sb, frags[sb].c_str());
    }

    if (opt.rfi && tot_cells)
        printf("RFI: flagged %.2f%% of (win,ant,pol) cells\n",
               100.0 * tot_flagged_cells / tot_cells);

    // ---- dedispersion shifts (integer native samples, ref = band top) ------
    std::vector<long long> shift(nchan_out, 0);
    if (opt.dm > 0.0) {
        const double ftop_ghz =
            freq_mhz(sb_lo == sb_hi ? sb_lo : 0, 0) * 1e-3;
        for (int ch = 0; ch < nchan_out; ch++) {
            const int sb = (sb_lo == sb_hi) ? sb_lo : ch / NCHAN_SB;
            const int c = (sb_lo == sb_hi) ? ch : ch % NCHAN_SB;
            const double f_ghz = freq_mhz(sb, c) * 1e-3;
            const double dt_ms =
                4.15 * opt.dm * (pow(f_ghz, -2.0) - pow(ftop_ghz, -2.0));
            shift[ch] = llround(dt_ms * 1e-3 / TSAMP_S);
        }
        printf("dedispersing to DM=%.3f (max shift %lld samples)\n", opt.dm,
               shift[nchan_out - 1]);
    }

    // ---- write SIGPROC ------------------------------------------------------
    FILE *fo = fopen(opt.out_fil.c_str(), "wb");
    if (!fo) {
        fprintf(stderr, "cannot open %s\n", opt.out_fil.c_str());
        return 1;
    }
    FilHeader hdr;
    hdr.source_name = opt.event.empty() ? "bbproc" : opt.event;
    hdr.nchans = nchan_out;
    hdr.fch1_mhz = freq_mhz(sb_lo == sb_hi ? sb_lo : 0, 0);
    hdr.foff_mhz = -DF_MHZ;
    hdr.tsamp_s = TSAMP_S * opt.tscrunch;
    hdr.tstart_mjd = tstart_mjd;
    hdr.telescope_id = opt.telescope_id;
    write_fil_header(fo, hdr);

    const long long nt_out = ntime / opt.tscrunch;
    std::vector<float> row(nchan_out);
    for (long long to = 0; to < nt_out; to++) {
        for (int ch = 0; ch < nchan_out; ch++) {
            float acc = 0.f;
            const long long tbase = to * opt.tscrunch + shift[ch];
            for (int k = 0; k < opt.tscrunch; k++) {
                const long long t = tbase + k;
                if (t >= 0 && t < ntime) acc += full[t * nchan_out + ch];
            }
            row[ch] = acc / opt.tscrunch;
        }
        fwrite(row.data(), 4, nchan_out, fo);
    }
    fclose(fo);
    printf("wrote %s: %lld samples x %d chans, tsamp %.3f us, tstart MJD "
           "%.9f\n", opt.out_fil.c_str(), nt_out, nchan_out,
           hdr.tsamp_s * 1e6, tstart_mjd);

    cudaFreeHost(h_block);
    cudaFree(d_raw);
    cudaFree(d_w);
    cudaFree(d_pow);
    if (d_autos) cudaFree(d_autos);
    if (d_mask) cudaFree(d_mask);
    if (d_wnorm) cudaFree(d_wnorm);
    return 0;
}

// ---------------------------------------------------------------------------
// Host: visibilities mode (legacy parity, single fragment)
// ---------------------------------------------------------------------------

struct VisOpts {
    std::string in_frag, out_vis, out_fil, cal_path, flag_path, delay_path;
    int tint = 8;
    bool averaging = false;
    int stokes = 0;
    float min_base = -1.f;
    long long npkts = -1, offpkts = 0;
    bool phase_only = false, swap_pol = false;
    int gpu = 0;
};

static int run_visibilities(const VisOpts &opt) {
    CUDA_CHECK(cudaSetDevice(opt.gpu));
    const long long fsz = file_size(opt.in_frag.c_str());
    if (fsz <= 0) {
        fprintf(stderr, "cannot open %s\n", opt.in_frag.c_str());
        return 1;
    }
    const long long tot_pkts = fsz / ((long long)NANT * NCHAN_SB * 4);
    const long long npkts =
        opt.npkts > 0 ? std::min(opt.npkts, tot_pkts - opt.offpkts)
                      : tot_pkts - opt.offpkts;
    printf("visibilities: %s, %lld packets (offset %lld), tint %d\n",
           opt.in_frag.c_str(), npkts, opt.offpkts, opt.tint);

    // cal weights in legacy per-coarse layout [ant][48][pol][ri]
    std::vector<float> wl((size_t)NANT * NCAL_COARSE * NPOL * 2);
    std::vector<float> antpos_e(NANT, 0.f);
    if (!opt.cal_path.empty()) {
        CalBlob cal;
        if (load_cal_blob(opt.cal_path.c_str(), &cal)) return 1;
        for (int a = 0; a < NANT; a++) {
            antpos_e[a] = cal.antpos_e[a];
            for (int c = 0; c < NCAL_COARSE; c++)
                for (int p = 0; p < NPOL; p++) {
                    const int cp = opt.swap_pol ? (1 - p) : p;
                    double gr = cal.gains[a][c][cp][0];
                    double gi = cal.gains[a][c][cp][1];
                    const double g2 = gr * gr + gi * gi;
                    if (g2 > 0.0) {
                        const double den = opt.phase_only ? sqrt(g2) : g2;
                        gr = gr / den;
                        gi = -gi / den;  // conj, matching legacy w = g*/|g|^2
                    }
                    wl[((size_t)a * NCAL_COARSE + c) * NPOL * 2 + p * 2] =
                        (float)gr;
                    wl[((size_t)a * NCAL_COARSE + c) * NPOL * 2 + p * 2 + 1] =
                        (float)gi;
                }
        }
    } else {
        for (int a = 0; a < NANT; a++)
            for (int c = 0; c < NCAL_COARSE; c++)
                for (int p = 0; p < NPOL; p++) {
                    wl[((size_t)a * NCAL_COARSE + c) * NPOL * 2 + p * 2] = 1.f;
                }
    }
    if (!opt.flag_path.empty()) {
        std::vector<int> fl;
        if (load_int_list(opt.flag_path.c_str(), &fl)) return 1;
        for (int a : fl)
            if (a >= 0 && a < NANT)
                memset(&wl[(size_t)a * NCAL_COARSE * NPOL * 2], 0,
                       NCAL_COARSE * NPOL * 2 * 4);
    }

    // baseline maps
    std::vector<int> a1(NBASE), a2(NBASE);
    int ctr = 0;
    for (int i = 0; i < NANT; i++)
        for (int j = i; j < NANT; j++) {
            a1[ctr] = i;
            a2[ctr] = j;
            ctr++;
        }
    std::vector<float> freqs(NCHAN_SB);
    const int sb = std::max(0, infer_sb_from_name(opt.in_frag));
    for (int c = 0; c < NCHAN_SB; c++) freqs[c] = (float)(freq_mhz(sb, c) * 1e6);

    std::vector<float> delays(NBASE, 0.f);
    const bool delaying = !opt.delay_path.empty();
    if (delaying) {
        FILE *fd = fopen(opt.delay_path.c_str(), "r");
        if (!fd) {
            fprintf(stderr, "cannot open %s\n", opt.delay_path.c_str());
            return 1;
        }
        for (int i = 0; i < NBASE; i++)
            if (fscanf(fd, "%f\n", &delays[i]) != 1) break;
        fclose(fd);
    }

    // device buffers
    unsigned char *h_block, *d_raw;
    float *d_prom, *d_corr, *d_final, *d_av, *d_fil, *d_wl, *d_freqs, *d_del,
        *d_ape;
    int *d_a1, *d_a2;
    CUDA_CHECK(cudaHostAlloc(&h_block, BLOCK_BYTES, cudaHostAllocDefault));
    CUDA_CHECK(cudaMalloc(&d_raw, BLOCK_BYTES));
    CUDA_CHECK(cudaMalloc(&d_prom, (size_t)NANT * NCHAN_SB * NPTR * 4));
    CUDA_CHECK(cudaMalloc(&d_corr, 2LL * NBASE * NCHAN_SB * 8 * 4));
    CUDA_CHECK(cudaMalloc(&d_final, (size_t)NBASE * NCHAN_SB * 8 * 4));
    CUDA_CHECK(cudaMalloc(&d_av, (size_t)NBASE * (NCHAN_SB / AV) * 4 * 4));
    CUDA_CHECK(cudaMalloc(&d_fil, NCHAN_SB * 4));
    CUDA_CHECK(cudaMalloc(&d_wl, wl.size() * 4));
    CUDA_CHECK(cudaMalloc(&d_freqs, NCHAN_SB * 4));
    CUDA_CHECK(cudaMalloc(&d_del, NBASE * 4));
    CUDA_CHECK(cudaMalloc(&d_ape, NANT * 4));
    CUDA_CHECK(cudaMalloc(&d_a1, NBASE * 4));
    CUDA_CHECK(cudaMalloc(&d_a2, NBASE * 4));
    CUDA_CHECK(cudaMemcpy(d_wl, wl.data(), wl.size() * 4,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_freqs, freqs.data(), NCHAN_SB * 4,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_del, delays.data(), NBASE * 4,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_ape, antpos_e.data(), NANT * 4,
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_a1, a1.data(), NBASE * 4, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_a2, a2.data(), NBASE * 4, cudaMemcpyHostToDevice));

    FILE *fin = fopen(opt.in_frag.c_str(), "rb");
    FILE *fout = opt.out_vis.empty() ? nullptr
                                     : fopen(opt.out_vis.c_str(), "wb");
    FILE *ffil = opt.out_fil.empty() ? nullptr
                                     : fopen(opt.out_fil.c_str(), "wb");
    if (!opt.out_vis.empty() && !fout) return 1;
    if (!opt.out_fil.empty() && !ffil) return 1;

    std::vector<float> outdata((size_t)NBASE * NCHAN_SB * 8);
    std::vector<float> filout(NCHAN_SB);
    const size_t pkt_bytes = (size_t)NANT * NCHAN_SB * 4;

    // seek to offset
    fseek(fin, opt.offpkts * pkt_bytes, SEEK_SET);

    auto flush_vis = [&](float *d_buf, float scf) {
        if (delaying)
            k_delay<<<NBASE * NCHAN_SB * 4 / 32, 32>>>(d_buf, d_freqs, d_del);
        if (fout) {
            if (opt.averaging) {
                k_fscrunch<<<NBASE * NCHAN_SB * 4 / AV / 32, 32>>>(d_buf, d_av);
                CUDA_CHECK(cudaMemcpy(outdata.data(), d_av,
                                      (size_t)NBASE * (NCHAN_SB / AV) * 4 * 4,
                                      cudaMemcpyDeviceToHost));
                fwrite(outdata.data(), 4, (size_t)NBASE * (NCHAN_SB / AV) * 4,
                       fout);
            } else {
                CUDA_CHECK(cudaMemcpy(outdata.data(), d_buf,
                                      outdata.size() * 4,
                                      cudaMemcpyDeviceToHost));
                fwrite(outdata.data(), 4, outdata.size(), fout);
            }
        }
        if (ffil) {
            k_reduce<<<NCHAN_SB, 1024>>>(d_buf, d_fil, scf, d_a1, d_a2,
                                         opt.stokes, d_ape, opt.min_base);
            CUDA_CHECK(cudaMemcpy(filout.data(), d_fil, NCHAN_SB * 4,
                                  cudaMemcpyDeviceToHost));
            fwrite(filout.data(), 4, NCHAN_SB, ffil);
        }
    };

    long long done = 0;
    int timi = 0;
    while (done < npkts) {
        const long long chunk_pkts =
            std::min((long long)PKT_PER_BLOCK, npkts - done);
        const size_t nread =
            fread(h_block, 1, chunk_pkts * pkt_bytes, fin);
        if (nread < pkt_bytes) break;
        const long long got = nread / pkt_bytes;
        CUDA_CHECK(cudaMemcpy(d_raw, h_block, got * pkt_bytes,
                              cudaMemcpyHostToDevice));
        for (long long p = 0; p < got; p++) {
            if (timi == 0) {
                k_zero<<<2 * NBASE * NCHAN_SB * 8 / 32, 32>>>(d_corr);
                k_zero<<<NBASE * NCHAN_SB * 8 / 32, 32>>>(d_final);
            }
            k_promote<<<NANT * NCHAN_SB * 4 / 32, 32>>>(
                d_raw + p * pkt_bytes, d_prom);
            k_correlate<<<NBASE * NCHAN_SB / 32 + 1, 32>>>(
                d_prom, d_corr, d_a1, d_a2, 1.f / opt.tint, d_wl);
            timi += 2;
            if (timi >= opt.tint) {
                if (opt.tint == 1) {
                    flush_vis(d_corr, 0.25f);
                    flush_vis(d_corr + (size_t)NBASE * NCHAN_SB * 8, 0.25f);
                } else {
                    k_adder<<<NBASE * NCHAN_SB * 8 / 32, 32>>>(d_corr, d_final);
                    flush_vis(d_final, 4.f);
                }
                timi = 0;
            }
        }
        done += got;
    }
    CUDA_CHECK(cudaDeviceSynchronize());
    fclose(fin);
    if (fout) fclose(fout);
    if (ffil) fclose(ffil);
    printf("visibilities done: %lld packets\n", done);
    return 0;
}

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

static void usage() {
    fprintf(stdout,
        "toolkit — DSA-110 M8 voltage-dump processing (dsa110-bbproc)\n"
        "\n"
        "Coherent filterbank mode:\n"
        "  -D <dir>          event directory (Level2/voltages or staging)\n"
        "  -E <event>        event name (fragments <event>_sbNN_data.out)\n"
        "  -i <file>         OR a single fragment (--sb to force subband)\n"
        "  -P <out.fil>      write full-band SIGPROC filterbank (float32)\n"
        "  --l <rad>         l offset from meridian phase centre [0]\n"
        "  --m <rad>         m offset [0]\n"
        "  -w <cal.dat>      beamformer_weights_*.dat cal blob\n"
        "  --phase-only      normalize cal gains to unit magnitude\n"
        "  --swap-pol        swap cal pol axis vs voltage pol axis\n"
        "  --core <file>     core antenna list [config/core_antennas.txt];\n"
        "                    'all' = every antenna\n"
        "  -f <file>         extra flagged antennas (voltage idx per line)\n"
        "  --tscrunch <n>    time integration [8] (262 us)\n"
        "  --dm <pc/cc>      integer-shift intra-channel dedispersion\n"
        "  --stokes <n>      0=I 1=Q 2=U 3=V [0]\n"
        "  --rfi             SK RFI flagging (realtime-equivalent stats)\n"
        "  --rfi-m <n>       SK window in native samples [256]\n"
        "  --mjd <mjd>       tstart override (else from manifest)\n"
        "  --telescope-id <n> SIGPROC telescope_id [0]\n"
        "\n"
        "Visibilities mode (legacy toolkit_dev parity, single fragment):\n"
        "  -i <file> -o <out.vis>  correlate all 4656 baselines\n"
        "  -t <n>            time integration in native samples [8]\n"
        "  -d <file>         per-baseline delays (ns, NBASE lines)\n"
        "  -a                8x frequency averaging (XX,YY only)\n"
        "  -p <out>          baseline-sum filterbank (legacy float stream)\n"
        "  -g <n>            its Stokes [0]\n"
        "  -v <m>            min E-W baseline length [none]\n"
        "  -s <n> -q <n>     packet count / offset\n"
        "\n"
        "Common:  --gpu <n> [0]   -h help\n");
}

int main(int argc, char *argv[]) {
    FilOpts fo;
    VisOpts vo;
    fo.core_path = "config/core_antennas.txt";
    bool want_fil = false, want_vis = false;

    enum { OPT_L = 1000, OPT_M, OPT_CORE, OPT_TSCR, OPT_DM, OPT_STOKES,
           OPT_RFI, OPT_RFIM, OPT_MJD, OPT_GPU, OPT_PHONLY, OPT_SWAP,
           OPT_SB, OPT_TELID };
    static struct option lopts[] = {
        {"l", required_argument, nullptr, OPT_L},
        {"m", required_argument, nullptr, OPT_M},
        {"core", required_argument, nullptr, OPT_CORE},
        {"tscrunch", required_argument, nullptr, OPT_TSCR},
        {"dm", required_argument, nullptr, OPT_DM},
        {"stokes", required_argument, nullptr, OPT_STOKES},
        {"rfi", no_argument, nullptr, OPT_RFI},
        {"rfi-m", required_argument, nullptr, OPT_RFIM},
        {"mjd", required_argument, nullptr, OPT_MJD},
        {"gpu", required_argument, nullptr, OPT_GPU},
        {"phase-only", no_argument, nullptr, OPT_PHONLY},
        {"swap-pol", no_argument, nullptr, OPT_SWAP},
        {"sb", required_argument, nullptr, OPT_SB},
        {"telescope-id", required_argument, nullptr, OPT_TELID},
        {nullptr, 0, nullptr, 0},
    };

    int c;
    while ((c = getopt_long(argc, argv, "D:E:i:P:w:f:o:t:d:ap:g:v:s:q:h",
                            lopts, nullptr)) != -1) {
        switch (c) {
            case 'D': fo.event_dir = optarg; break;
            case 'E': fo.event = optarg; break;
            case 'i': fo.single_frag = optarg; vo.in_frag = optarg; break;
            case 'P': fo.out_fil = optarg; want_fil = true; break;
            case 'w': fo.cal_path = optarg; vo.cal_path = optarg; break;
            case 'f': fo.flag_path = optarg; vo.flag_path = optarg; break;
            case 'o': vo.out_vis = optarg; want_vis = true; break;
            case 't': vo.tint = atoi(optarg); break;
            case 'd': vo.delay_path = optarg; break;
            case 'a': vo.averaging = true; break;
            case 'p': vo.out_fil = optarg; want_vis = true; break;
            case 'g': vo.stokes = atoi(optarg); break;
            case 'v': vo.min_base = atof(optarg); break;
            case 's': vo.npkts = atoll(optarg); break;
            case 'q': vo.offpkts = atoll(optarg); break;
            case OPT_L: fo.l = atof(optarg); break;
            case OPT_M: fo.m = atof(optarg); break;
            case OPT_CORE: fo.core_path = optarg; break;
            case OPT_TSCR: fo.tscrunch = atoi(optarg); break;
            case OPT_DM: fo.dm = atof(optarg); break;
            case OPT_STOKES: fo.stokes = atoi(optarg); break;
            case OPT_RFI: fo.rfi = true; break;
            case OPT_RFIM: fo.rfi_m = atoi(optarg); break;
            case OPT_MJD: fo.mjd_override = atof(optarg); break;
            case OPT_GPU: fo.gpu = atoi(optarg); vo.gpu = atoi(optarg); break;
            case OPT_PHONLY: fo.phase_only = true; vo.phase_only = true; break;
            case OPT_SWAP: fo.swap_pol = true; vo.swap_pol = true; break;
            case OPT_SB: fo.single_sb = atoi(optarg); break;
            case OPT_TELID: fo.telescope_id = atoi(optarg); break;
            case 'h': usage(); return 0;
            default: usage(); return 1;
        }
    }

    if (want_fil && want_vis) {
        fprintf(stderr, "pick one mode: -P (filterbank) or -o/-p "
                "(visibilities)\n");
        return 1;
    }
    if (want_fil) {
        if (fo.event_dir.empty() == fo.single_frag.empty()) {
            fprintf(stderr, "filterbank mode needs -D/-E or -i\n");
            return 1;
        }
        if (!fo.event_dir.empty() && fo.event.empty()) {
            fprintf(stderr, "-D needs -E <event>\n");
            return 1;
        }
        if (fo.tscrunch < 1 || (T_PER_BLOCK % fo.tscrunch)) {
            fprintf(stderr, "--tscrunch must divide %d\n", T_PER_BLOCK);
            return 1;
        }
        return run_filterbank(fo);
    }
    if (want_vis) {
        if (vo.in_frag.empty()) {
            fprintf(stderr, "visibilities mode needs -i <fragment>\n");
            return 1;
        }
        return run_visibilities(vo);
    }
    usage();
    return 1;
}
