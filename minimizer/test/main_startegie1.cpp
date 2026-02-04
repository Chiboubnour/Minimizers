
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

#include "code.hpp"   // t_minimizer_sequence_19_16

// ================= CONFIG =================
#define CHUNK_BASES   (1 << 20)        // 1 MB
#define TARGET_MB  8
#define TARGET_BASES  (TARGET_MB * 1024 * 1024ULL)

#define MAX_HASHES    (CHUNK_BASES)
#define FREQ_HZ       100e6
#define SEQ_FILE_PATH "0:/data"

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

// ================= PACK ASCII → uint64 =================
void pack_bases_ascii(
    const uint8_t* seq,
    uint64_t* packed,
    uint32_t n_bases
) {
    uint32_t words = (n_bases + 7) / 8;
    for (uint32_t i = 0; i < words; i++) packed[i] = 0;

    for (uint32_t i = 0; i < n_bases; i++) {
        uint32_t w = i >> 3;
        uint32_t b = i & 7;
        packed[w] |= ((uint64_t)seq[i]) << (8 * b);
    }
}

// ================= MAIN =================
int main() {

    xil_printf("\n============================================\n\r");
    xil_printf(" 🚀 SD → CPU vs FPGA Minimizer (FIXED SIZE)\n\r");
    xil_printf("============================================\n\r");

    xil_printf("Target dataset : %llu bases (%.2f MB)\n\r",
               TARGET_BASES,
               (double)TARGET_BASES / (1024.0 * 1024.0));

    // -------- SD INIT --------
    if (init_sd_card() != XST_SUCCESS) {
        xil_printf("❌ SD init failed\n\r");
        return -1;
    }

    FIL file;
    if (f_open(&file, SEQ_FILE_PATH, FA_READ) != FR_OK) {
        xil_printf("❌ Cannot open FASTA\n\r");
        return -1;
    }

    // -------- BUFFERS --------
    uint8_t*  seq_ascii  = (uint8_t*)  malloc(CHUNK_BASES);
    uint64_t* seq_packed = (uint64_t*) malloc(((CHUNK_BASES + 7) / 8) * sizeof(uint64_t));
    uint64_t* hash_cpu   = (uint64_t*) malloc(MAX_HASHES * sizeof(uint64_t));
    uint64_t* hash_hw    = (uint64_t*) malloc(MAX_HASHES * sizeof(uint64_t));

    if (!seq_ascii || !seq_packed || !hash_cpu || !hash_hw) {
        xil_printf("❌ Allocation failed\n\r");
        return -1;
    }

    // -------- INIT IP --------
    XMinimizer accel;
    auto* cfg = XMinimizer_LookupConfig(XPAR_MINIMIZER_0_DEVICE_ID);
    XMinimizer_CfgInitialize(&accel, cfg);

    // -------- STATS --------
    bool eof = false;
    uint64_t remaining_bases = TARGET_BASES;

    uint64_t total_bases_cpu = 0;
    uint64_t total_bases_hw  = 0;
    uint64_t total_min_cpu   = 0;
    uint64_t total_min_hw    = 0;

    double total_time_cpu = 0.0;
    double total_time_hw  = 0.0;

    // -------- LOOP --------
    while (!eof && remaining_bases > 0) {

        uint32_t max_to_read =
            (remaining_bases < CHUNK_BASES)
            ? (uint32_t)remaining_bases
            : CHUNK_BASES;

        uint32_t n_bases = read_sequence_chunk(
            file, seq_ascii, max_to_read, eof
        );

        if (n_bases == 0)
            break;

        remaining_bases -= n_bases;

        pack_bases_ascii(seq_ascii, seq_packed, n_bases);

        // ===== CPU =====
        XTime t1, t2;
        XTime_GetTime(&t1);

        int n_cpu = t_minimizer_sequence_19_16(
            (const char*)seq_ascii,
            n_bases,
            hash_cpu
        );

        XTime_GetTime(&t2);
        double t_cpu = (double)(t2 - t1) / FREQ_HZ;

        total_bases_cpu += n_bases;
        total_min_cpu   += n_cpu;
        total_time_cpu  += t_cpu;

        // ===== FPGA =====
        uint32_t words = (n_bases + 7) / 8;

        Xil_DCacheFlushRange((UINTPTR)seq_packed, words * sizeof(uint64_t));
        Xil_DCacheFlushRange((UINTPTR)hash_hw, MAX_HASHES * sizeof(uint64_t));

        XMinimizer_Set_packed_sequence(&accel, (UINTPTR)seq_packed);
         XMinimizer_Set_out_hashes(&accel, (UINTPTR)hash_hw);
         XMinimizer_Set_n_bases(&accel, n_bases);
        XTime t3, t4;
        XTime_GetTime(&t3);

        XMinimizer_Start(&accel);
        while (!XMinimizer_IsDone(&accel));

        XTime_GetTime(&t4);
        double t_hw = (double)(t4 - t3) / FREQ_HZ;

        Xil_DCacheInvalidateRange((UINTPTR)hash_hw, MAX_HASHES * sizeof(uint64_t));
        uint32_t n_hw = XMinimizer_Get_n_minimizers(&accel);

        total_bases_hw += n_bases;
        total_min_hw   += n_hw;
        total_time_hw  += t_hw;
    }

    // -------- CLEANUP --------
    f_close(&file);
    free(seq_ascii);
    free(seq_packed);
    free(hash_cpu);
    free(hash_hw);

    // -------- RESULTS --------
    xil_printf("\n=====================================\n\r");
    xil_printf(" 📊 GLOBAL STATISTICS (%.2f MB)\n\r",
               (double)total_bases_cpu / (1024.0 * 1024.0));
    xil_printf("=====================================\n\r");

    xil_printf("[CPU]\n\r");
    xil_printf("Bases traitées     : %llu\n\r", total_bases_cpu);
    xil_printf("Minimizers générés : %llu\n\r", total_min_cpu);
    printf("Temps total CPU    : %.6f s\n\r", total_time_cpu);

    xil_printf("\n[FPGA]\n\r");
    xil_printf("Bases traitées     : %llu\n\r", total_bases_hw);
    xil_printf("Minimizers générés : %llu\n\r", total_min_hw);
    printf("Temps total FPGA   : %.6f s\n\r", total_time_hw);

    if (total_time_hw > 0.0) {
        printf("\n[SPEEDUP] CPU / FPGA = x%.2f\n\r",
                   total_time_cpu / total_time_hw);
    }

    double thr_cpu = (double)total_bases_cpu / total_time_cpu / 1e6;
    double thr_hw  = (double)total_bases_hw  / total_time_hw  / 1e6;

    xil_printf("\n[THROUGHPUT]\n\r");
    printf("CPU  : %.2f Mbases/s\n\r", thr_cpu);
    printf("FPGA : %.2f Mbases/s\n\r", thr_hw);

    xil_printf("\n=== DONE ✅ ===\n\r");
    return 0;
}

