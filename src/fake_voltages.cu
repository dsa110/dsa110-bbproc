// -*- c++ -*-
// fake_voltages — synthesize an M8-format DSA-110 voltage-dump event
// (16 subband fragments + manifests) containing thermal int4 noise plus
// a dispersed pulse at a given (l, m) offset, mirroring the dsa110-rt
// injection convention (dsart/inject/online.py):
//
//   v_ant(t, ch) = noise + A * box(t - t_arr(ch); W)
//                          * e^{+2*pi*i * nu * (E*l + N*m)/c} [* g_ant]
//
//   t_arr(ch) = t0 + 4.15 ms * DM * (f_GHz^-2 - ftop_GHz^-2)
//
// The pulse amplitude is deterministic per (ant, ch) — the same model
// the realtime OnlineInjector uses — so a coherent beamformer pointed
// at (l, m) recovers amplitude^2 * N_ant^2 while an offset pointing
// decorrelates. With -w the signal is multiplied by the cal-blob gains
// so the toolkit's cal application round-trips too.
//
// Noise is a deterministic counter-hash Gaussian (Box-Muller over a
// splitmix64 stream), quantized to int4 — no curand dependency.

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <getopt.h>

#include <string>
#include <vector>

#include "bbproc.h"

using namespace bbproc;

// ---- deterministic per-cell RNG -------------------------------------------
__device__ __host__ inline unsigned long long splitmix64(unsigned long long x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// two iid N(0,1) from one 64-bit counter
__device__ inline void gauss2(unsigned long long ctr, float *g1, float *g2) {
    const unsigned long long a = splitmix64(ctr);
    const unsigned long long b = splitmix64(ctr ^ 0xDEADBEEFCAFEF00DULL);
    // u in (0,1]
    const float u1 = ((a >> 11) + 1) * (1.0f / 9007199254740992.0f);
    const float u2 = (b >> 11) * (1.0f / 9007199254740992.0f);
    const float r = sqrtf(-2.0f * logf(u1));
    const float th = 2.0f * (float)M_PI * u2;
    *g1 = r * cosf(th);
    *g2 = r * sinf(th);
}

__device__ inline unsigned char pack_int4(float re, float im) {
    int r = __float2int_rn(re);
    int i = __float2int_rn(im);
    r = max(-8, min(7, r));
    i = max(-8, min(7, i));
    return (unsigned char)((r & 0x0F) | ((i & 0x0F) << 4));
}

// One thread per (pkt, ant); loops channels. Writes 4 bytes per (ch).
// sig[ant][ch][2]: complex signal factor A * phasor (* gain), same for
// both pols. tarr[ch]: pulse arrival (native samples, global). width in
// native samples. blk0_t: global native-sample index of this block's start.
__global__ void k_fake_block(unsigned char *__restrict__ out,
                             const float *__restrict__ sig,
                             const long long *__restrict__ tarr, int width,
                             long long blk0_t, float noise_sigma,
                             unsigned long long seed) {
    const int pkt = blockIdx.x;
    const int ant = threadIdx.x + blockIdx.y * blockDim.x;
    if (ant >= NANT) return;
    for (int ch = 0; ch < NCHAN_SB; ch++) {
        const float sr = sig[((size_t)ant * NCHAN_SB + ch) * 2];
        const float si = sig[((size_t)ant * NCHAN_SB + ch) * 2 + 1];
        const size_t base =
            (((size_t)pkt * NANT + ant) * NCHAN_SB + ch) * NT_PER_PKT * NPOL;
        for (int tt = 0; tt < NT_PER_PKT; tt++) {
            const long long t = blk0_t + (long long)pkt * NT_PER_PKT + tt;
            const bool on = (t >= tarr[ch]) && (t < tarr[ch] + width);
            for (int p = 0; p < NPOL; p++) {
                const unsigned long long ctr =
                    seed ^ ((unsigned long long)t << 22) ^
                    ((unsigned long long)(ant * NCHAN_SB + ch) << 2) ^
                    (unsigned long long)p;
                float g1, g2;
                gauss2(ctr, &g1, &g2);
                float re = noise_sigma * g1, im = noise_sigma * g2;
                if (on) {
                    re += sr;
                    im += si;
                }
                out[base + (size_t)tt * NPOL + p] = pack_int4(re, im);
            }
        }
    }
}

static void usage() {
    fprintf(stdout,
        "fake_voltages — synthesize an M8-format voltage event\n"
        "  -O <dir>        output directory [.]\n"
        "  -E <event>      event name [fake0000test]\n"
        "  --nblocks <n>   blocks per fragment [4] (23 = production)\n"
        "  --l <rad>       source l offset [0]\n"
        "  --m <rad>       source m offset [0]\n"
        "  --dm <pc/cc>    dispersion measure [100]\n"
        "  --width <n>     pulse width, native samples [16]\n"
        "  --amp <v>       per-ant per-pol signal amplitude, int4 units [1.0]\n"
        "  --t0 <s>        pulse arrival at band top, s from file start [0.5]\n"
        "  --noise <v>     noise sigma in int4 units [1.33]\n"
        "  -w <cal.dat>    cal blob: antpos for the phasor + gains folded\n"
        "                  into the signal (omit: antpos=0, unit gains)\n"
        "  --mjd <mjd>     manifest mjd_target [60000.0]\n"
        "  --dec-deg <deg> imprint the F21 dec fringe (meridian source at\n"
        "                  this pointing dec); recover with toolkit --dec-deg\n"
        "  --sb <a>-<b>    subband range [0-15]\n"
        "  --seed <n>      RNG seed [12345]\n"
        "  --gpu <n>       [0]\n");
}

int main(int argc, char *argv[]) {
    std::string outdir = ".", event = "fake0000test", cal_path;
    int nblocks = 4, width = 16, gpu = 0, sb_lo = 0, sb_hi = NSB - 1;
    double l = 0.0, m = 0.0, dm = 100.0, amp = 1.0, t0 = 0.5, noise = 1.33;
    double mjd = 60000.0, dec_deg = NAN;
    unsigned long long seed = 12345;

    enum { OPT_NB = 1000, OPT_L, OPT_M, OPT_DM, OPT_W, OPT_AMP, OPT_T0,
           OPT_NOISE, OPT_MJD, OPT_SB, OPT_SEED, OPT_GPU, OPT_DEC };
    static struct option lopts[] = {
        {"nblocks", required_argument, nullptr, OPT_NB},
        {"l", required_argument, nullptr, OPT_L},
        {"m", required_argument, nullptr, OPT_M},
        {"dm", required_argument, nullptr, OPT_DM},
        {"width", required_argument, nullptr, OPT_W},
        {"amp", required_argument, nullptr, OPT_AMP},
        {"t0", required_argument, nullptr, OPT_T0},
        {"noise", required_argument, nullptr, OPT_NOISE},
        {"mjd", required_argument, nullptr, OPT_MJD},
        {"dec-deg", required_argument, nullptr, OPT_DEC},
        {"sb", required_argument, nullptr, OPT_SB},
        {"seed", required_argument, nullptr, OPT_SEED},
        {"gpu", required_argument, nullptr, OPT_GPU},
        {nullptr, 0, nullptr, 0},
    };
    int c;
    while ((c = getopt_long(argc, argv, "O:E:w:h", lopts, nullptr)) != -1) {
        switch (c) {
            case 'O': outdir = optarg; break;
            case 'E': event = optarg; break;
            case 'w': cal_path = optarg; break;
            case OPT_NB: nblocks = atoi(optarg); break;
            case OPT_L: l = atof(optarg); break;
            case OPT_M: m = atof(optarg); break;
            case OPT_DM: dm = atof(optarg); break;
            case OPT_W: width = atoi(optarg); break;
            case OPT_AMP: amp = atof(optarg); break;
            case OPT_T0: t0 = atof(optarg); break;
            case OPT_NOISE: noise = atof(optarg); break;
            case OPT_MJD: mjd = atof(optarg); break;
            case OPT_DEC: dec_deg = atof(optarg); break;
            case OPT_SB: sscanf(optarg, "%d-%d", &sb_lo, &sb_hi); break;
            case OPT_SEED: seed = strtoull(optarg, nullptr, 10); break;
            case OPT_GPU: gpu = atoi(optarg); break;
            case 'h': usage(); return 0;
            default: usage(); return 1;
        }
    }
    if (event.size() > 16) {
        fprintf(stderr, "event name > 16 chars breaks the C2 wire "
                "convention; refusing\n");
        return 1;
    }

    CUDA_CHECK(cudaSetDevice(gpu));

    CalBlob cal;
    bool have_cal = false;
    if (!cal_path.empty()) {
        if (load_cal_blob(cal_path.c_str(), &cal)) return 1;
        have_cal = true;
        printf("using antpos + gains from %s\n", cal_path.c_str());
    }

    const double ftop_ghz = freq_mhz(0, 0) * 1e-3;
    constexpr double LAT_OVRO_RAD = 0.6498558936875687;
    const double sin_dec_lat =
        std::isnan(dec_deg) ? 0.0
                            : sin(dec_deg * M_PI / 180.0 - LAT_OVRO_RAD);
    if (!std::isnan(dec_deg))
        printf("dec fringe: dec=%.4f deg sin(dec-lat)=%.6f\n",
               dec_deg, sin_dec_lat);
    printf("event %s: %d blocks/frag, DM %.2f, width %d, amp %.2f, "
           "(l,m)=(%.6g,%.6g), t0 %.3f s, sb %d..%d\n",
           event.c_str(), nblocks, dm, width, amp, l, m, t0, sb_lo, sb_hi);

    unsigned char *d_out, *h_block;
    float *d_sig;
    long long *d_tarr;
    CUDA_CHECK(cudaMalloc(&d_out, BLOCK_BYTES));
    CUDA_CHECK(cudaHostAlloc(&h_block, BLOCK_BYTES, cudaHostAllocDefault));
    CUDA_CHECK(cudaMalloc(&d_sig, (size_t)NANT * NCHAN_SB * 2 * 4));
    CUDA_CHECK(cudaMalloc(&d_tarr, NCHAN_SB * 8));

    std::vector<float> sig((size_t)NANT * NCHAN_SB * 2);
    std::vector<long long> tarr(NCHAN_SB);

    for (int sb = sb_lo; sb <= sb_hi; sb++) {
        // per-(ant, ch) complex signal factor
        for (int a = 0; a < NANT; a++) {
            const double E = have_cal ? cal.antpos_e[a] : 0.0;
            const double N = have_cal ? cal.antpos_n[a] : 0.0;
            for (int ch = 0; ch < NCHAN_SB; ch++) {
                const double f_hz = freq_mhz(sb, ch) * 1e6;
                // A raw-voltage source at (l, m) off the dec-stopped
                // meridian carries: the (l,m) phasor (+ sign, SNAP
                // convention), PLUS the dec fringe
                // e^{-2 pi i f sin(dec-lat) N/c} (dsart cal_loader),
                // PLUS conj(cal weight) — the toolkit multiplies the
                // blob in unconjugated, so v must carry its conjugate
                // for the product to come out flat. Getting all three
                // right here is what lets tests/roundtrip.sh catch
                // convention bugs (260715twmx incident).
                // sign pinned EMPIRICALLY by 260715twmx: the toolkit
                // weight e^{-2 pi i f (El+Nm+sinD*N)/c} recovers the
                // real FRB, so raw voltages carry the + sign of all
                // three terms.
                const double geo =
                    2.0 * M_PI * f_hz *
                    (E * l + N * m + sin_dec_lat * N) / CVAC;
                double sr = amp * cos(geo), si = amp * sin(geo);  // +i geo
                if (have_cal) {
                    // conj(pol-B weight) into both pols (test
                    // simplification; per-pol only matters for QUV)
                    const double gr = cal.gains[a][ch / FINE_PER_COARSE][0][0];
                    const double gi = -cal.gains[a][ch / FINE_PER_COARSE][0][1];
                    const double g2 = gr * gr + gi * gi;
                    if (g2 > 0.0) {
                        const double tr = (sr * gr - si * gi) / sqrt(g2);
                        si = (sr * gi + si * gr) / sqrt(g2);
                        sr = tr;
                    }
                }
                sig[((size_t)a * NCHAN_SB + ch) * 2] = (float)sr;
                sig[((size_t)a * NCHAN_SB + ch) * 2 + 1] = (float)si;
            }
        }
        // arrival times
        for (int ch = 0; ch < NCHAN_SB; ch++) {
            const double f_ghz = freq_mhz(sb, ch) * 1e-3;
            const double dt_s = 4.15e-3 * dm *
                (pow(f_ghz, -2.0) - pow(ftop_ghz, -2.0));
            tarr[ch] = llround((t0 + dt_s) / TSAMP_S);
        }
        CUDA_CHECK(cudaMemcpy(d_sig, sig.data(), sig.size() * 4,
                              cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_tarr, tarr.data(), NCHAN_SB * 8,
                              cudaMemcpyHostToDevice));

        char sbtag[8];
        snprintf(sbtag, sizeof sbtag, "_sb%02d", sb);
        const std::string dpath =
            outdir + "/" + event + sbtag + "_data.out";
        FILE *f = fopen(dpath.c_str(), "wb");
        if (!f) {
            fprintf(stderr, "cannot open %s\n", dpath.c_str());
            return 1;
        }
        for (int b = 0; b < nblocks; b++) {
            dim3 grid(PKT_PER_BLOCK, (NANT + 95) / 96);
            k_fake_block<<<grid, 96>>>(d_out, d_sig, d_tarr, width,
                                       (long long)b * T_PER_BLOCK,
                                       (float)noise, seed + sb * 1000003ULL);
            CUDA_CHECK(cudaMemcpy(h_block, d_out, BLOCK_BYTES,
                                  cudaMemcpyDeviceToHost));
            fwrite(h_block, 1, BLOCK_BYTES, f);
        }
        fclose(f);

        // manifest (fields the toolkit + C3 collection care about)
        const std::string mpath = outdir + "/" + event + sbtag + ".json";
        FILE *mf = fopen(mpath.c_str(), "w");
        if (mf) {
            fprintf(mf,
                "{\n"
                "  \"event_name\": \"%s\",\n"
                "  \"cn_id\": %d,\n"
                "  \"chgroup\": %d,\n"
                "  \"subband\": \"sb%02d\",\n"
                "  \"synthetic\": true,\n"
                "  \"target_block_n\": %d,\n"
                "  \"block_n_first\": 0,\n"
                "  \"block_n_last\": %d,\n"
                "  \"n_blocks_written\": %d,\n"
                "  \"n_blocks_dropped\": 0,\n"
                "  \"bytes_per_block\": %lld,\n"
                "  \"total_bytes\": %lld,\n"
                "  \"mjd_target\": %.9f,\n"
                "  \"fake\": {\"l\": %.9g, \"m\": %.9g, \"dm\": %.4f,\n"
                "           \"width\": %d, \"amp\": %.4f, \"t0_s\": %.6f,\n"
                "           \"noise_sigma\": %.4f, \"seed\": %llu}\n"
                "}\n",
                event.c_str(), sb, sb, sb, 0, nblocks - 1, nblocks,
                BLOCK_BYTES, (long long)nblocks * BLOCK_BYTES,
                mjd, l, m, dm, width, amp, t0, noise, seed);
            fclose(mf);
        }
        printf("wrote %s (%d blocks)\n", dpath.c_str(), nblocks);
    }
    cudaFree(d_out);
    cudaFree(d_sig);
    cudaFree(d_tarr);
    cudaFreeHost(h_block);
    return 0;
}
