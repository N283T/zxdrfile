// C reference benchmark for xdrfile (mdtraj).
// Compares throughput against Zig implementation.
//
// Build: make -C benchmarks/c_reference
// Run:   ./benchmarks/c_reference/bench
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "xdrfile_trr.h"

static long get_file_size(const char *path) {
    struct stat st;
    if (stat(path, &st) != 0) return -1;
    return st.st_size;
}

static double time_ms(struct timespec *start, struct timespec *end) {
    return (end->tv_sec - start->tv_sec) * 1000.0 +
           (end->tv_nsec - start->tv_nsec) / 1000000.0;
}

static void bench_xtc(const char *path, const char *name) {
    int natoms = 0;
    if (read_xtc_natoms((char *)path, &natoms) != 0) return;

    long file_size = get_file_size(path);
    if (file_size < 0) return;

    XDRFILE *xd = xdrfile_open(path, "r");
    if (!xd) return;

    rvec *x = (rvec *)malloc(natoms * sizeof(rvec));
    matrix box;
    int step;
    float time_val, prec;
    int n_frames = 0;

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    while (read_xtc(xd, natoms, &step, &time_val, box, x, &prec) == 0) {
        n_frames++;
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);

    double elapsed = time_ms(&t0, &t1);
    double size_mb = file_size / (1024.0 * 1024.0);
    double throughput = size_mb / (elapsed / 1000.0);
    double fps = n_frames / (elapsed / 1000.0);

    printf("  %-30s %6d atoms  %5d frames  %7.1f MB  %8.1f ms  %7.1f MB/s  %8.0f fps\n",
           name, natoms, n_frames, size_mb, elapsed, throughput, fps);

    free(x);
    xdrfile_close(xd);
}

static void bench_trr(const char *path, const char *name) {
    int natoms = 0;
    if (read_trr_natoms((char *)path, &natoms) != 0) return;

    long file_size = get_file_size(path);
    if (file_size < 0) return;

    XDRFILE *xd = xdrfile_open(path, "r");
    if (!xd) return;

    rvec *x = (rvec *)malloc(natoms * sizeof(rvec));
    rvec *v = (rvec *)malloc(natoms * sizeof(rvec));
    rvec *f = (rvec *)malloc(natoms * sizeof(rvec));
    matrix box;
    int step;
    float time_val, lambda;
    int n_frames = 0;

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    while (read_trr(xd, natoms, &step, &time_val, &lambda, box, x, v, f) == 0) {
        n_frames++;
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);

    double elapsed = time_ms(&t0, &t1);
    double size_mb = file_size / (1024.0 * 1024.0);
    double throughput = size_mb / (elapsed / 1000.0);
    double fps = n_frames / (elapsed / 1000.0);

    printf("  %-30s %6d atoms  %5d frames  %7.1f MB  %8.1f ms  %7.1f MB/s  %8.0f fps\n",
           name, natoms, n_frames, size_mb, elapsed, throughput, fps);

    free(x);
    free(v);
    free(f);
    xdrfile_close(xd);
}

int main(void) {
    printf("\n=== C xdrfile benchmark (mdtraj) ===\n\n");

    printf("XTC:\n");
    bench_xtc("test_data/1l2y.xtc", "1l2y (small)");
    bench_xtc("benchmarks/md_data/3tvj_I_R1.xtc", "3tvj_I (531 atoms)");
    bench_xtc("benchmarks/md_data/5wvo_C_R1.xtc", "5wvo_C (3858 atoms)");
    bench_xtc("benchmarks/md_data/6sup_A_R1.xtc", "6sup_A (33377 atoms)");

    printf("\nTRR:\n");
    bench_trr("test_data/frame0.trr", "frame0 (small)");
    bench_trr("benchmarks/md_data/3tvj_I_R1.trr", "3tvj_I (531 atoms)");
    bench_trr("benchmarks/md_data/5wvo_C_R1.trr", "5wvo_C (3858 atoms)");
    bench_trr("benchmarks/md_data/6sup_A_R1.trr", "6sup_A (33377 atoms)");

    printf("\n");
    return 0;
}
