/*
#include "xparameters.h"
#include "xil_printf.h"
#include "xil_cache.h"
#include "xminimizer.h"
#include "xminimizer_store.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>


static const char TEST_FASTA[] =
    ">Seq1\n"
    "ACCTCCTTTATGGAGGACTACACCTATCACTTCGAGAAAGGTAATGACCTGGTGCTCGGCTCCCATATGCTGGAAGTCTGCCCGTCGATCGCCGCAGAAG\n"
    ">Seq2\n"
    "CCAGATAGTCAAAAGAGACATGCCGTTGTACTATCTGAACTCTTCGTGGCTTAAGAATCGTAGAACCATCGGCTGATACATTCTGAACAAGAATATCTAC\n";

#define RAW_SIZE ((uint32_t)(sizeof(TEST_FASTA) - 1))   // -1 : exclut le '\0' final

struct hash128_t {
    uint64_t lo;
    uint64_t hi;
} __attribute__((packed, aligned(16)));

#define BURST_LEN 512
#define MAX_HASH_SLOTS (2 * BURST_LEN)

int main() {

    xil_printf("\n============================================\n\r");
    xil_printf("Test piste 2 : (PARSER+COMPUTE fusionnes) -> STORE (2 IP) sur FASTA fixe\n\r");
    xil_printf("Build : %s %s\n\r", __DATE__, __TIME__);
    xil_printf("============================================\n\r");

    xil_printf("FASTA brut (%lu octets) :\n%s\n\r", (unsigned long)RAW_SIZE, TEST_FASTA);

    // -------- BUFFERS --------
    // raw_fasta : texte FASTA brut (entete + sequence), tel quel -- lu
    //             directement par le kernel fusionne, aucun buffer
    //             intermediaire de bases filtrees necessaire.
    // hash_hw   : sortie finale de store (paires 128 bits)

    uint8_t*   raw_fasta = (uint8_t*)   malloc(RAW_SIZE);
    hash128_t* hash_hw   = (hash128_t*) malloc(MAX_HASH_SLOTS * sizeof(hash128_t));

    if (!raw_fasta || !hash_hw) {
        xil_printf("Allocation failed\n\r");
        return -1;
    }

    memcpy(raw_fasta, TEST_FASTA, RAW_SIZE);


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

    // -------- EXECUTION --------

    Xil_DCacheFlushRange((UINTPTR)raw_fasta, RAW_SIZE);
    Xil_DCacheFlushRange((UINTPTR)hash_hw, MAX_HASH_SLOTS * sizeof(hash128_t));

    XMinimizer_Set_raw_ptr(&accel_minimizer, (UINTPTR)raw_fasta);
    XMinimizer_Set_n_bytes(&accel_minimizer, RAW_SIZE);

    XMinimizer_store_Set_tab_hash(&accel_store, (UINTPTR)hash_hw);
    // store ne connait pas a l'avance le nombre de bases FILTREES (le
    // filtrage se fait desormais a l'interieur du kernel fusionne, jamais
    // expose au CPU) -- on lui passe le nombre d'octets BRUTS comme
    // majorant sur (n_bases+GROUP_W-1)/GROUP_W (borne de securite pour le
    // dimensionnement de la boucle d'ecriture en rafale) ; l'ecriture
    // s'arrete de toute facon des que le flux signale sa vraie fin
    // (comportement deja inchange dans thr_write_burst).
    XMinimizer_store_Set_n_bases(&accel_store, RAW_SIZE);

    XMinimizer_store_Start(&accel_store);
    XMinimizer_Start(&accel_minimizer);

    while (!XMinimizer_store_IsDone(&accel_store));
    while (!XMinimizer_IsDone(&accel_minimizer));

    uint64_t n_hw = XMinimizer_store_Get_nMinizrs(&accel_store);

    // -------- AFFICHAGE DES RESULTATS --------

    Xil_DCacheInvalidateRange((UINTPTR)hash_hw, MAX_HASH_SLOTS * sizeof(hash128_t));

    xil_printf("\n=====================================\n\r");
    printf("Minimizers trouves : %llu\n", (unsigned long long)n_hw);
    xil_printf("=====================================\n\r");

    for (uint64_t i = 0; i < n_hw; i++) {
        uint64_t group_idx = i / 2;
        uint64_t val = (i % 2 == 0) ? hash_hw[group_idx].lo : hash_hw[group_idx].hi;
        printf("Minimizer [%4llu] : 0x%016llx\n", (unsigned long long)i, (unsigned long long)val);
    }

    xil_printf("\n=====================================\n\r");

    free(raw_fasta);
    free(hash_hw);

    return 0;
}
*/
