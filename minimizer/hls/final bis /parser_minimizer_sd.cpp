#include "xparameters.h"
#include "xil_printf.h"
#include "xil_cache.h"
#include "xminimizer.h"
#include "xminimizer_store.h"
#include "xtime_l.h"
#include "xsdps.h"
#include "ff.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define CHUNK_SIZE_BYTES (8ULL * 1024 * 1024)   // taille de chunk FIXE : 8 Mo de FASTA brut par chunk

#define TOTAL_SIZE_MB    32   // <-- taille TOTALE a extraire du fichier SD : modifier cette valeur pour tester d'autres volumes
#define TOTAL_SIZE_BYTES (TOTAL_SIZE_MB * 1024ULL * 1024ULL)

#define SEQ_FILE_PATH "0:/data"

struct hash128_t {
    uint64_t lo;
    uint64_t hi;
} __attribute__((packed, aligned(16)));

// Majorant large : au pire tout le chunk brut n'est que des bases valides.
#define MAX_HASH_SLOTS ((CHUNK_SIZE_BYTES / 2) + 1)

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

int main() {

    xil_printf("\n============================================\n\r");
    xil_printf("SD -> FPGA Minimizer Benchmark piste 2 (2 IP : minimizer[parser+compute] + store)\n\r");
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

    uint8_t*   raw_chunk = (uint8_t*)   malloc(CHUNK_SIZE_BYTES);
    hash128_t* hash_hw   = (hash128_t*) malloc(MAX_HASH_SLOTS * sizeof(hash128_t));

    if (!raw_chunk || !hash_hw) {
        xil_printf("Allocation failed\n\r");
        return -1;
    }


    XMinimizer accel_minimizer;
    XMinimizer_Config *cfg_minimizer =
        XMinimizer_LookupConfig(XPAR_MINIMIZER_0_DEVICE_ID);

    if (!cfg_minimizer || XMinimizer_CfgInitialize(&accel_minimizer, cfg_minimizer) != XST_SUCCESS) {
        xil_printf("IP minimizer (parser+compute) init failed\n\r");
        return -1;
    }

    XMinimizer_store accel_store;
    XMinimizer_store_Config *cfg_store =
        XMinimizer_store_LookupConfig(XPAR_MINIMIZER_STORE_0_DEVICE_ID);

    if (!cfg_store || XMinimizer_store_CfgInitialize(&accel_store, cfg_store) != XST_SUCCESS) {
        xil_printf("IP store init failed\n\r");
        return -1;
    }

    xil_printf("IP minimizer (parser+compute) + store initialisees\n\r");

    // -------- STATS --------

    uint64_t total_raw_bytes   = 0;
    uint64_t total_minimizers  = 0;
    double   total_kernel_time = 0.0;
    uint64_t chunk_idx         = 0;

    UINT br;


    const uint64_t n_chunks = (TOTAL_SIZE_BYTES + CHUNK_SIZE_BYTES - 1) / CHUNK_SIZE_BYTES;

    printf("Taille totale a traiter     : %llu octets (%d Mo)\n",
           (unsigned long long)TOTAL_SIZE_BYTES, TOTAL_SIZE_MB);
    printf("Taille de chunk (fixe)      : %llu octets (8 Mo)\n",
           (unsigned long long)CHUNK_SIZE_BYTES);
    printf("Nombre de chunks a traiter  : %llu\n", (unsigned long long)n_chunks);

    // -------- LOOP (sequentielle, sans chevauchement de chunks) --------

    for (uint64_t c = 0; c < n_chunks; c++) {

        const uint64_t remaining = TOTAL_SIZE_BYTES - total_raw_bytes;
        const uint64_t to_read   = (remaining < CHUNK_SIZE_BYTES) ? remaining : CHUNK_SIZE_BYTES;

        if (f_read(&file, raw_chunk, to_read, &br) != FR_OK) {
            xil_printf("Erreur de lecture SD\n\r");
            break;
        }
        if (br == 0) break;   // fin de fichier (SD plus courte que TOTAL_SIZE_MB)

        total_raw_bytes += br;

        Xil_DCacheFlushRange((UINTPTR)raw_chunk, br);
        Xil_DCacheFlushRange((UINTPTR)hash_hw, MAX_HASH_SLOTS * sizeof(hash128_t));

        XMinimizer_Set_raw_ptr(&accel_minimizer, (UINTPTR)raw_chunk);
        XMinimizer_Set_n_bytes(&accel_minimizer, br);

        XMinimizer_store_Set_tab_hash(&accel_store, (UINTPTR)hash_hw);
        XMinimizer_store_Set_n_bases(&accel_store, br);

        XTime k1, k2;
        XTime_GetTime(&k1);

        XMinimizer_store_Start(&accel_store);
        XMinimizer_Start(&accel_minimizer);

        while (!XMinimizer_store_IsDone(&accel_store));

        XTime_GetTime(&k2);

        while (!XMinimizer_IsDone(&accel_minimizer));

        total_kernel_time += (double)(k2 - k1) / freq_hz;

        uint64_t n_hw = XMinimizer_store_Get_nMinizrs(&accel_store);
        total_minimizers += n_hw;

        if (chunk_idx == 0) {

            Xil_DCacheInvalidateRange((UINTPTR)hash_hw, MAX_HASH_SLOTS * sizeof(hash128_t));

            uint64_t n_show = (n_hw < 5) ? n_hw : 5;
            printf("5 premiers hash trouves (chunk 0) :\n");
            for (uint64_t i = 0; i < n_show; i++) {
                uint64_t group_idx = i / 2;
                uint64_t val = (i % 2 == 0) ? hash_hw[group_idx].lo : hash_hw[group_idx].hi;
                printf("  Minimizer [%2llu] : 0x%016llx\n", (unsigned long long)i, (unsigned long long)val);
            }
        }

        chunk_idx++;
    }

    // -------- CLEANUP --------

    f_close(&file);

    free(raw_chunk);
    free(hash_hw);

    // =================================================
    // THROUGHPUT
    // =================================================

    double thr_kernel = (double)total_raw_bytes / total_kernel_time / 1e6;
    double ns_per_bp  = total_kernel_time * 1e9 / (double)total_raw_bytes;

    xil_printf("\n=====================================\n\r");
    xil_printf("GLOBAL RESULTS piste 2 (2 IP : minimizer[parser+compute] + store)\n\r");
    xil_printf("=====================================\n\r");
    printf("Chunks traites              : %llu\n", (unsigned long long)chunk_idx);
    printf("Octets bruts lus            : %llu\n", (unsigned long long)total_raw_bytes);
    printf("Total minimizers (HW)       : %llu\n", (unsigned long long)total_minimizers);
    printf("Temps minimizer+store       : %.6f s\n", total_kernel_time);
    printf("Debit FPGA (approx, brut)   : %.2f Moctets-bruts/s\n", thr_kernel);
    printf("Temps par base (approx)     : %.4f ns/bp\n", ns_per_bp);
    xil_printf("\n=====================================\n\r");

    return 0;
}
