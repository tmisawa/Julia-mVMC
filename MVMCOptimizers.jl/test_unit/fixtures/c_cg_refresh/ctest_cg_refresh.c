/* Optional fixed-input oracle: verbatim C 1.4.0 CG Mains/operator, serial only.
 * See ctest_cg_refresh.md. No Cargo/build-script dependency.
 * Usage: probe <input> <output> [diag]   (diag: DiagScale-preconditioned Main) */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int MPI_Comm;
#define MPI_DOUBLE 0
#ifdef MVMC_SRCG_REAL
#define USE_IMAG 0
#else
#define USE_IMAG 1
#endif
static int NVMCSample, NSROptCGMaxIter;
static double DSROptCGTol = 0.0, DSROptStaDel = 1e-5, Wc;
static int MPI_Comm_rank(MPI_Comm c, int *r) { (void)c; *r = 0; return 0; }
static int MPI_Comm_size(MPI_Comm c, int *s) { (void)c; *s = 1; return 0; }
static int MPI_Bcast(double *x, int n, int t, int r, MPI_Comm c) {
    (void)x; (void)n; (void)t; (void)r; (void)c; return 0;
}
static int MPI_Barrier(MPI_Comm c) { (void)c; return 0; }
static void SafeMpiAllReduce(double *a, double *b, int n, MPI_Comm c) {
    (void)c; memcpy(b, a, (size_t)n * sizeof(*b));
}
#define StartTimer(n) ((void)(n))
#define StopTimer(n) ((void)(n))
extern void dgemv_(const char *, const int *, const int *, const double *,
                   const double *, const int *, const double *, const int *,
                   const double *, double *, const int *);
#define M_DGEMV dgemv_
/* GNU inline linkage only; the extracted sequential arithmetic is unchanged. */
#define inline static
#include "ctest_cg_dot_upstream.inc"
#undef inline
int fn_operate_by_S(int, double *, double *, double *, MPI_Comm);
#include "ctest_cg_main_upstream.inc"

static int read_bits(FILE *f, double *v, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        uint64_t bits;
        if (fscanf(f, "%" SCNx64, &bits) != 1) return 0;
        memcpy(v + i, &bits, sizeof(bits));
    }
    return 1;
}
static void write_bits(FILE *f, const double *v, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        uint64_t bits;
        memcpy(&bits, v + i, sizeof(bits));
        fprintf(f, "%s%016" PRIx64, i ? " " : "", bits);
    }
    fputc('\n', f);
}
int main(int argc, char **argv) {
    if (argc != 3 && argc != 4) return 2;
    int precon = argc == 4 && strcmp(argv[3], "diag") == 0;
    if (argc == 4 && !precon) return 2;
    FILE *input = fopen(argv[1], "r");
    if (!input) return 3;
    int ch;
    while ((ch = fgetc(input)) == '#') {
        while ((ch = fgetc(input)) != '\n' && ch != EOF) {}
    }
    if (ch != EOF) ungetc(ch, input);
    int n, complex;
    if (fscanf(input, "%d %d %d", &n, &NVMCSample, &complex) != 3 ||
        n < 1 || n > 10000 || NVMCSample < 1 || NVMCSample > 10000 ||
        complex != USE_IMAG) return 4;
    Wc = (double)NVMCSample;
    size_t count = (size_t)n * 8 + (size_t)NVMCSample * (USE_IMAG + 1) * (n + 1);
    double *base = calloc(count, sizeof(*base));
    double *trial = calloc(count, sizeof(*trial));
    if (!base || !trial) return 5;
    double *mean = base + 3*n, *diag = base + 2*n, *real = base + 4*n;
    double *imag = real + (size_t)n*NVMCSample, *gradient = base + n;
    if (!read_bits(input, mean, n) || !read_bits(input, diag, n) ||
        !read_bits(input, real, (size_t)n*NVMCSample) ||
        !read_bits(input, imag, (size_t)USE_IMAG*n*NVMCSample) ||
        !read_bits(input, gradient, n)) return 5;
    fclose(input); /* Only independently archived input records are consumed. */
    FILE *out = fopen(argv[2], "w");
    if (!out) return 6;
    fprintf(out, "# C 1.4.0 %s/operator verbatim; serial LP64 OpenBLAS; tol=0 shift=1e-5\n",
            precon ? "StochasticOptPreconCG_DiagScale_Main" : "StochasticOptCG_Main");
    fprintf(out, "%d %d %d\n", n, NVMCSample, complex);
    write_bits(out, mean, n); write_bits(out, diag, n);
    write_bits(out, real, (size_t)n*NVMCSample);
    write_bits(out, imag, (size_t)USE_IMAG*n*NVMCSample);
    write_bits(out, gradient, n);
    double *local = imag + (size_t)USE_IMAG*n*NVMCSample + (USE_IMAG + 1)*NVMCSample;
    fn_operate_by_S(n, gradient, base, base, 0);
    write_bits(out, base, n);
    memset(base, 0, (size_t)n*sizeof(*base));
    size_t d_offset = (size_t)(local - base) + 2*n;
    for (int limit = 1; limit <= 41; ++limit) {
        memcpy(trial, base, count*sizeof(*base));
        NSROptCGMaxIter = limit;
        int iterations = precon ? fn_StochasticOptPreconCG_DiagScale_Main(n, trial, 0)
                              : fn_StochasticOptCG_Main(n, trial, 0);
        fprintf(out, "%d %d ", limit, iterations);
        write_bits(out, trial, n);
        write_bits(out, trial + d_offset + n, n);
        write_bits(out, trial + d_offset, n);
    }
    free(base); free(trial);
    return fclose(out) == 0 ? 0 : 7;
}
