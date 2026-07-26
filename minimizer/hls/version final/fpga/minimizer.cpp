
#include "xparameters.h"
#include "xil_printf.h"
#include "xil_cache.h"
#include "xminimizer.h"
#include "xtime_l.h"
#include "xsdps.h"
#include "ff.h"

#include <cstdint>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>


#define CHUNK_SIZE_BYTES (1048576ULL)

#define TOTAL_SIZE_MB 8
#define TOTAL_SIZE_BYTES (TOTAL_SIZE_MB * 1024ULL * 1024ULL)

#define TARGET_BASES TOTAL_SIZE_BYTES
#define CHUNK_BASES  CHUNK_SIZE_BYTES

#define SEQ_FILE_PATH "0:/data"

// tab_hash : paires de 128 bits (2 x uint64_t) desormais.
struct hash128_t {
    uint64_t lo;
    uint64_t hi;
} __attribute__((packed, aligned(16)));

#define MAX_HASH_SLOTS ((CHUNK_BASES / 2) + 1)


// ================= FATFS =================
static FATFS fatfs;

int init_sd_card() {
    XSdPs_Config *cfg;
    XSdPs sd;
    int status;

    xil_printf("[SD] Initialisation...\n\r");

    cfg = XSdPs_LookupConfig(XPAR_XSDPS_0_DEVICE_ID);
    if (!cfg) return XST_FAILURE;

    status = XSdPs_CfgInitialize(&sd, cfg, cfg->BaseAddress);
    if (status != XST_SUCCESS) return XST_FAILURE;

    if (f_mount(&fatfs, "0:/", 1) != FR_OK) return XST_FAILURE;

    xil_printf("[SD] OK\n\r");
    return XST_SUCCESS;
}

// ================= READ FASTA =================
uint32_t read_sequence_chunk(
    FIL& file,
    uint8_t* buffer,
    uint32_t max_bases,
    bool& eof
) {
    UINT br;
    char c;
    uint32_t n = 0;
    eof = false;

    while (n < max_bases) {
        if (f_read(&file, &c, 1, &br) != FR_OK || br == 0) {
            eof = true;
            break;
        }

        if (c == '>') {
            do {
                if (f_read(&file, &c, 1, &br) != FR_OK || br == 0) {
                    eof = true;
                    break;
                }
            } while (c != '\n');
            continue;
        }

        if (c == '\n' || c == '\r') continue;

        c = toupper(c);
        if (c=='A' || c=='C' || c=='G' || c=='T')
            buffer[n++] = (uint8_t)c;
    }
    return n;
}

int main() {

    xil_printf("\n============================================\n\r");
    xil_printf("SD -> FPGA Minimizer Benchmark (kernel unique, debit seul)\n\r");
    xil_printf("Build : %s %s\n\r", __DATE__, __TIME__);
    xil_printf("============================================\n\r");

    const double freq_hz = (double)COUNTS_PER_SECOND;

    // -------- SD INIT --------
    if (init_sd_card() != XST_SUCCESS) {
        xil_printf("SD init failed\n\r");
        return -1;
    }

    FIL file;
    if (f_open(&file, SEQ_FILE_PATH, FA_READ) != FR_OK) {
        xil_printf("Cannot open FASTA\n\r");
        return -1;
    }

    const uint64_t total_chunks = (TOTAL_SIZE_BYTES + CHUNK_SIZE_BYTES - 1) / CHUNK_SIZE_BYTES;

    xil_printf("\n[INFO] Dataset size : %llu MB\n\r", (unsigned long long)TOTAL_SIZE_MB);
    xil_printf("[INFO] Chunk size   : %llu MB\n\r", (unsigned long long)(CHUNK_SIZE_BYTES/(1024*1024)));
    xil_printf("[INFO] Total chunks : %llu\n\r", (unsigned long long)total_chunks);

    // -------- BUFFERS --------

    uint8_t*   seq_ascii = (uint8_t*)   malloc(CHUNK_SIZE_BYTES);
    hash128_t* hash_hw   = (hash128_t*) malloc(MAX_HASH_SLOTS * sizeof(hash128_t));

    if (!seq_ascii || !hash_hw) {
        xil_printf("Allocation failed\n\r");
        return -1;
    }

    // -------- INIT IP --------

    XMinimizer accel;
    XMinimizer_Config *cfg = XMinimizer_LookupConfig(XPAR_MINIMIZER_0_DEVICE_ID);

    if (!cfg || XMinimizer_CfgInitialize(&accel, cfg) != XST_SUCCESS) {
        xil_printf("IP init failed\n\r");
        return -1;
    }

    xil_printf("IP initialisee\n\r");

    // -------- STATS --------

    bool eof = false;
    uint64_t remaining_bases = TARGET_BASES;
    uint64_t total_bases = 0;
    double fpga_kernel_time = 0.0;

    uint64_t total_min_hw  = 0;
    uint64_t chunk_idx = 0;

    // -------- LOOP --------

    while (!eof && remaining_bases > 0) {

        uint32_t max_to_read = (remaining_bases < CHUNK_BASES) ? (uint32_t)remaining_bases : CHUNK_BASES;

        uint32_t n_bases = read_sequence_chunk(file, seq_ascii, max_to_read, eof);

        if (n_bases == 0)
            break;

        remaining_bases -= n_bases;
        total_bases += n_bases;

        // =================================================
        // FPGA EXECUTION (kernel unique)
        // =================================================

        Xil_DCacheFlushRange((UINTPTR)seq_ascii, n_bases);
        Xil_DCacheFlushRange((UINTPTR)hash_hw, MAX_HASH_SLOTS * sizeof(hash128_t));

        XMinimizer_Set_packed_sequence(&accel, (UINTPTR)seq_ascii);
        XMinimizer_Set_tab_hash(&accel, (UINTPTR)hash_hw);
        XMinimizer_Set_n_bases(&accel, n_bases);

        XTime k1, k2;

        XTime_GetTime(&k1);

        XMinimizer_Start(&accel);

        while (!XMinimizer_IsDone(&accel));

        XTime_GetTime(&k2);

        double kernel_time = (double)(k2 - k1) / freq_hz;
        fpga_kernel_time += kernel_time;

        uint64_t n_hw = XMinimizer_Get_nMinizrs(&accel);
        total_min_hw += n_hw;

        chunk_idx++;
    }

    // -------- CLEANUP --------

    f_close(&file);

    free(seq_ascii);
    free(hash_hw);

    // =================================================
    // THROUGHPUT
    // =================================================

    double thr_kernel = (double)total_bases / fpga_kernel_time / 1e6;

    xil_printf("\n=====================================\n\r");
    xil_printf("GLOBAL RESULTS (kernel unique)\n\r");
    xil_printf("=====================================\n\r");

    printf("\n--- EXECUTION TIME ---\n");
    printf("FPGA kernel time (total) : %.6f s\n", fpga_kernel_time);

    printf("\n--- THROUGHPUT ---\n");
    printf("FPGA kernel throughput   : %.2f Mbases/s\n", thr_kernel);

    printf("\n--- RESULTATS ---\n");
    printf("Total bases traitees      : %llu\n", (unsigned long long)total_bases);
    printf("Total minimizers (HW)     : %llu\n", (unsigned long long)total_min_hw);
    printf("Chunks traites            : %llu\n", (unsigned long long)chunk_idx);

    xil_printf("\n=====================================\n\r");

    return 0;
}
