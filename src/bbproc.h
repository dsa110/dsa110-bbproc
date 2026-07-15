// -*- c++ -*-
// dsa110-bbproc: shared definitions for offline baseband processing of
// DSA-110 M8 voltage dumps (dsa110-rt `voltage_retention` fragments).
//
// Voltage fragment format (per corr node / subband, NO header):
//   <event>_sb<NN>_data.out = N_block x 288 MiB fada blocks, each
//   [2048 packets, 96 ants, 384 chans, 2 times, 2 pols] int4-complex
//   bytes: real = low nibble, imag = high nibble (sign-extended int4).
//   One packet = 2 native time samples of 32.768 us.
//   (dsa110-rt: dsart/dump/voltage_ring.py, slow_corr_kernel.unpack_int4_split)
//
// Frequency plan (dsa110-rt dsart/common/constants.py freq_GHz):
//   f(g, c) = 1530.0 MHz - (1024 + g*384 + c) * (250/8192) MHz
//   Subband g=0 channel 0 (band top) = 1498.75 MHz, descending.
//
// Cal blob "beamformer_weights_*.dat" (74,496 bytes, float32 LE;
// dsa110-rt dsart/cal/bf_weights.py):
//   antpos_e[96] (m), antpos_n[96] (m),
//   gains[96 ant][48 coarse ch][2 pol][re, im]   (pol order [B, A])
//   Each coarse channel applies to 8 adjacent fine channels.
//
// (l, m) phasor convention (dsa110-rt dsart/inject/online.py, "DSA-110
// SNAP convention", POSITIVE sign): a source offset (l, m) rad from the
// meridian phase centre appears in calibrated voltages with
//   e^{+2*pi*i * nu * (E*l + N*m + U*n) / c},   n = sqrt(1-l^2-m^2) - 1.
// Beamforming applies the conjugate. U ~ 0 for the (approximately
// planar) core; outriggers are excluded from the coherent sum anyway.

#ifndef BBPROC_H
#define BBPROC_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>

namespace bbproc {

// ---- array / format constants (pinned to dsa110-rt) -----------------------
constexpr int NANT = 96;
constexpr int NCHAN_SB = 384;          // channels per subband (chgroup)
constexpr int NSB = 16;                // subbands = corr nodes
constexpr int NCHAN_FULL = NSB * NCHAN_SB;   // 6144
constexpr int NPOL = 2;
constexpr int NT_PER_PKT = 2;          // native samples per packet
constexpr int PKT_PER_BLOCK = 2048;
constexpr long long BLOCK_BYTES =
    (long long)PKT_PER_BLOCK * NANT * NCHAN_SB * NT_PER_PKT * NPOL; // 301,989,888
constexpr int T_PER_BLOCK = PKT_PER_BLOCK * NT_PER_PKT;             // 4096
constexpr double TSAMP_S = 32.768e-6;  // native sample
constexpr double F0_CONF_MHZ = 1530.0;
constexpr double DF_MHZ = 250.0 / 8192.0;              // 0.030517578125
constexpr int CH0_SYS = 1024;                          // system chan of sb0 ch0
constexpr double FTOP_MHZ = F0_CONF_MHZ - CH0_SYS * DF_MHZ;  // 1498.75
constexpr double CVAC = 299792458.0;

// cal blob
constexpr int NCAL_COARSE = 48;                        // per subband
constexpr int FINE_PER_COARSE = NCHAN_SB / NCAL_COARSE;  // 8
constexpr long long CAL_BLOB_BYTES =
    4LL * (2 * NANT + NANT * NCAL_COARSE * NPOL * 2);  // 74,496

inline double freq_mhz(int sb, int ch) {
    return F0_CONF_MHZ - (CH0_SYS + sb * NCHAN_SB + ch) * DF_MHZ;
}

// ---- cal blob --------------------------------------------------------------
struct CalBlob {
    float antpos_e[NANT];
    float antpos_n[NANT];
    // gains[ant][coarse][pol] complex (re, im), pol order [B, A]
    float gains[NANT][NCAL_COARSE][NPOL][2];
};

// Load a beamformer_weights_*.dat blob. Returns 0 on success.
inline int load_cal_blob(const char *path, CalBlob *out) {
    FILE *f = fopen(path, "rb");
    if (!f) { fprintf(stderr, "cal blob: cannot open %s\n", path); return 1; }
    fseek(f, 0, SEEK_END);
    long long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (sz != CAL_BLOB_BYTES) {
        fprintf(stderr, "cal blob: %s is %lld bytes, expected %lld\n",
                path, sz, CAL_BLOB_BYTES);
        fclose(f);
        return 1;
    }
    size_t n = 0;
    n += fread(out->antpos_e, sizeof(float), NANT, f);
    n += fread(out->antpos_n, sizeof(float), NANT, f);
    n += fread(out->gains, sizeof(float), NANT * NCAL_COARSE * NPOL * 2, f);
    fclose(f);
    return n == (size_t)(2 * NANT + NANT * NCAL_COARSE * NPOL * 2) ? 0 : 1;
}

// ---- small text loaders ----------------------------------------------------
// One integer per line; '#' comments and blank lines ignored.
inline int load_int_list(const char *path, std::vector<int> *out) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return 1; }
    char line[256];
    while (fgets(line, sizeof line, f)) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\0') continue;
        out->push_back(atoi(p));
    }
    fclose(f);
    return 0;
}

// ---- minimal manifest scraping ---------------------------------------------
// The staging manifest (<event>_sb<NN>.json) is flat json; we only need a
// couple of numeric fields, so scan for '"key": <number>' without a json lib.
inline bool json_get_double(const std::string &text, const char *key,
                            double *out) {
    std::string pat = std::string("\"") + key + "\"";
    size_t p = text.find(pat);
    if (p == std::string::npos) return false;
    p = text.find(':', p + pat.size());
    if (p == std::string::npos) return false;
    *out = strtod(text.c_str() + p + 1, nullptr);
    return true;
}

inline bool read_text_file(const char *path, std::string *out) {
    FILE *f = fopen(path, "r");
    if (!f) return false;
    char buf[4096];
    size_t n;
    while ((n = fread(buf, 1, sizeof buf, f)) > 0) out->append(buf, n);
    fclose(f);
    return true;
}

// ---- SIGPROC filterbank header --------------------------------------------
inline void sig_put_string(FILE *f, const char *s) {
    int32_t n = (int32_t)strlen(s);
    fwrite(&n, sizeof n, 1, f);
    fwrite(s, 1, n, f);
}
inline void sig_put_int(FILE *f, const char *k, int32_t v) {
    sig_put_string(f, k); fwrite(&v, sizeof v, 1, f);
}
inline void sig_put_double(FILE *f, const char *k, double v) {
    sig_put_string(f, k); fwrite(&v, sizeof v, 1, f);
}
inline void sig_put_str(FILE *f, const char *k, const char *v) {
    sig_put_string(f, k); sig_put_string(f, v);
}

struct FilHeader {
    std::string source_name = "unknown";
    double fch1_mhz = FTOP_MHZ;   // centre freq of first (highest) channel
    double foff_mhz = -DF_MHZ;
    int nchans = NCHAN_FULL;
    double tsamp_s = TSAMP_S;
    double tstart_mjd = 0.0;
    int nbits = 32;
    int nifs = 1;
    int telescope_id = 0;         // no official DSA-110 id; override via CLI
    int machine_id = 0;
    double src_raj = 0.0;         // hhmmss.s convention (optional)
    double src_dej = 0.0;         // ddmmss.s convention (optional)
};

inline void write_fil_header(FILE *f, const FilHeader &h) {
    sig_put_string(f, "HEADER_START");
    sig_put_str(f, "source_name", h.source_name.c_str());
    sig_put_int(f, "telescope_id", h.telescope_id);
    sig_put_int(f, "machine_id", h.machine_id);
    sig_put_int(f, "data_type", 1);            // filterbank
    sig_put_double(f, "fch1", h.fch1_mhz);
    sig_put_double(f, "foff", h.foff_mhz);
    sig_put_int(f, "nchans", h.nchans);
    sig_put_int(f, "nbits", h.nbits);
    sig_put_double(f, "tstart", h.tstart_mjd);
    sig_put_double(f, "tsamp", h.tsamp_s);
    sig_put_int(f, "nifs", h.nifs);
    if (h.src_raj != 0.0) sig_put_double(f, "src_raj", h.src_raj);
    if (h.src_dej != 0.0) sig_put_double(f, "src_dej", h.src_dej);
    sig_put_string(f, "HEADER_END");
}

// ---- fragment paths ---------------------------------------------------------
inline std::string frag_path(const std::string &dir, const std::string &event,
                             int sb) {
    char buf[32];
    snprintf(buf, sizeof buf, "_sb%02d_data.out", sb);
    return dir + "/" + event + buf;
}
inline std::string manifest_path(const std::string &dir,
                                 const std::string &event, int sb) {
    char buf[32];
    snprintf(buf, sizeof buf, "_sb%02d.json", sb);
    return dir + "/" + event + buf;
}

inline long long file_size(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    fseek(f, 0, SEEK_END);
    long long s = ftell(f);
    fclose(f);
    return s;
}

// ---- CUDA error check -------------------------------------------------------
#define CUDA_CHECK(call)                                                     \
    do {                                                                     \
        cudaError_t _e = (call);                                             \
        if (_e != cudaSuccess) {                                             \
            fprintf(stderr, "CUDA error %s at %s:%d: %s\n", #call,           \
                    __FILE__, __LINE__, cudaGetErrorString(_e));             \
            exit(1);                                                         \
        }                                                                    \
    } while (0)

}  // namespace bbproc

#endif  // BBPROC_H
