#include "xparameters.h"
#include "xil_printf.h"
#include "xil_cache.h"
#include "xminimizer_compute.h"
#include "xminimizer_store.h"
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

// tab_hash : paires de 128 bits (2 x uint64_t) -- format inchange, ecrit par
// le kernel "store" (thr_write_burst), identique au kernel unique.
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
    xil_printf("SD -> FPGA Minimizer Benchmark (2 IP : compute + store, debit seul)\n\r");
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


    // -------- BUFFERS --------

    uint8_t*   seq_ascii = (uint8_t*)   malloc(CHUNK_SIZE_BYTES);
    hash128_t* hash_hw   = (hash128_t*) malloc(MAX_HASH_SLOTS * sizeof(hash128_t));

    if (!seq_ascii || !hash_hw) {
        xil_printf("Allocation failed\n\r");
        return -1;
    }

    // -------- INIT IP (2 instances distinctes) --------
    //
    // NOTE : XPAR_MINIMIZER_COMPUTE_0_DEVICE_ID / XPAR_MINIMIZER_STORE_0_DEVICE_ID
    // et les noms XMinimizer_compute_* / XMinimizer_store_* supposent que les
    // instances Vivado s'appellent "minimizer_compute_0" / "minimizer_store_0"
    // (noms par defaut donnes par create_bd_cell). A ajuster si tu les as
    // renommees -- verifier les macros exactes dans xparameters.h une fois la
    // plateforme Vitis regeneree depuis le nouveau XSA (2 IP + FIFO).

    XMinimizer_compute accel_compute;
    XMinimizer_compute_Config *cfg_compute =
        XMinimizer_compute_LookupConfig(XPAR_MINIMIZER_COMPUTE_0_DEVICE_ID);

    if (!cfg_compute || XMinimizer_compute_CfgInitialize(&accel_compute, cfg_compute) != XST_SUCCESS) {
        xil_printf("IP compute init failed\n\r");
        return -1;
    }

    XMinimizer_store accel_store;
    XMinimizer_store_Config *cfg_store =
        XMinimizer_store_LookupConfig(XPAR_MINIMIZER_STORE_0_DEVICE_ID);

    if (!cfg_store || XMinimizer_store_CfgInitialize(&accel_store, cfg_store) != XST_SUCCESS) {
        xil_printf("IP store init failed\n\r");
        return -1;
    }

    xil_printf("IP compute + store initialisees\n\r");

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

        if (chunk_idx == 0) {

            printf("Premiere base lue : %c\n", (char)seq_ascii[0]);
        }

        Xil_DCacheFlushRange((UINTPTR)seq_ascii, n_bases);
        Xil_DCacheFlushRange((UINTPTR)hash_hw, MAX_HASH_SLOTS * sizeof(hash128_t));
        XMinimizer_compute_Set_packed_sequence(&accel_compute, (UINTPTR)seq_ascii);
        XMinimizer_compute_Set_n_bases(&accel_compute, n_bases);

        XMinimizer_store_Set_tab_hash(&accel_store, (UINTPTR)hash_hw);
        XMinimizer_store_Set_n_bases(&accel_store, n_bases);

        XTime k1, k2;

        XTime_GetTime(&k1);

        XMinimizer_store_Start(&accel_store);
        XMinimizer_compute_Start(&accel_compute);


        while (!XMinimizer_store_IsDone(&accel_store));

        XTime_GetTime(&k2);

        while (!XMinimizer_compute_IsDone(&accel_compute));

        double kernel_time = (double)(k2 - k1) / freq_hz;
        fpga_kernel_time += kernel_time;

        uint64_t n_hw = XMinimizer_store_Get_nMinizrs(&accel_store);
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
    xil_printf("GLOBAL RESULTS (2 IP : compute + store)\n\r");
    xil_printf("=====================================\n\r");
    xil_printf("\n[INFO] Dataset size : %llu MB\n\r", (unsigned long long)TOTAL_SIZE_MB);
    printf("Total bases traitees   : %llu\n", (unsigned long long)total_bases);
    printf("Total minimizers (HW)  : %llu\n", (unsigned long long)total_min_hw);
    printf("FPGA kernel time       : %.6f s\n", fpga_kernel_time);
    printf("FPGA kernel throughput : %.2f Mbases/s\n", thr_kernel);
    xil_printf("\n=====================================\n\r");

    return 0;
}
