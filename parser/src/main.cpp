#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cinttypes>
#include <ap_int.h>
#include "functions.hpp"

#if 0
int main(int argc, char * argv[])
{
    printf("=====================================\n");
    printf(" Minimizer test   : %s\n", __FILE__);
    printf(" Compilation date : %s %s\n", __DATE__, __TIME__);
    printf("=====================================\n");

    printf("\n✅ DONE\n");
    return 0;
}
#endif

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <immintrin.h>

#define BUF_SIZE 65536

/* ASCII -> 2 bits
   A=0 C=1 G=2 T=3
   autres = 255 */
static uint8_t lut[256];

void init_lut()
{
    for(int i=0;i<256;i++)
        lut[i]=255;

    lut['A']=0; lut['a']=0;
    lut['C']=1; lut['c']=1;
    lut['G']=2; lut['g']=2;
    lut['T']=3; lut['t']=3;
}

/* Parser FASTA */
void parse_fasta(FILE *in, FILE *out)
{
    uint8_t buffer[] = ">seq1\nACGTACGTACGT\nNNNN\n>seq2\nTTTTCCCCAAAAGGGG";

    int in_header = 0;

    uint8_t packed = 0;
    int shift = 6;

    size_t n;

    while((n=fread(buffer,1,BUF_SIZE,in))>0)
    {
        size_t i=0;

        /* --- SIMD 32 bytes --- */
        for(; i+32<=n; i+=32)
        {
            __m256i v =
                _mm256_loadu_si256(
                    (__m256i*)(buffer+i));

            uint8_t tmp[32];
            _mm256_storeu_si256(
                (__m256i*)tmp,v);

            for(int j=0;j<32;j++)
            {
                uint8_t c = tmp[j];

                if(c=='>')
                {
                    in_header=1;
                    continue;
                }

                if(c=='\n')
                {
                    in_header=0;
                    continue;
                }

                if(in_header)
                    continue;

                uint8_t code = lut[c];
                if(code==255)
                    continue;

                packed |= (code<<shift);
                shift -=2;

                if(shift<0)
                {
                    fwrite(&packed,1,1,out);
                    packed=0;
                    shift=6;
                }
            }
        }

        /* --- reste scalaire --- */
        for(;i<n;i++)
        {
            uint8_t c = buffer[i];

            if(c=='>')
            {
                in_header=1;
                continue;
            }

            if(c=='\n')
            {
                in_header=0;
                continue;
            }

            if(in_header)
                continue;

            uint8_t code = lut[c];
            if(code==255)
                continue;

            packed |= (code<<shift);
            shift-=2;

            if(shift<0)
            {
                fwrite(&packed,1,1,out);
                packed=0;
                shift=6;
            }
        }
    }

    /* flush final */
    if(shift!=6)
        fwrite(&packed,1,1,out);
}

int main(int argc,char**argv)
{
    if(argc!=3)
    {
        printf("usage: %s input.fasta output.bin\n",
               argv[0]);
        return 1;
    }

    FILE *in=fopen(argv[1],"rb");
    if(!in){perror("input");return 1;}

    FILE *out=fopen(argv[2],"wb");
    if(!out){perror("output");return 1;}

    init_lut();
    parse_fasta(in,out);

    fclose(in);
    fclose(out);

    printf("Conversion OK\n");
    return 0;
}
