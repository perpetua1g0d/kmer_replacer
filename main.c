#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>



#define mask_size 1 << 16
#define kmer_max 0b11111111
#define extended_bits 0b11111100


uint8_t dna_encode(const char c)
{
    if (c == 'A')
        return 0b00;
    if (c == 'C')
        return 0b01;
    if (c == 'G')
        return 0b10;
    if (c == 'T')
        return 0b11;

    return -1;
}

char dna_decode(const uint16_t x)
{
    if (x == 0b00)
        return 'A';
    if (x == 0b01)
        return 'C';
    if (x == 0b10)
        return 'G';
    if (x == 0b11)
        return 'T';

    return -1;
}

void dna_kmer_decode(uint16_t dna_kmer, char *s)
{
    for (int i = 0; i < 8; ++i)
    {
        s[i] = dna_decode(dna_kmer & 0b11);
        dna_kmer >>= 2;
    }
    s[8] = '\0';
}

uint8_t iupac_encode(const char c)
{
    if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
        return dna_encode(c);

    if (c == 'R')
        return dna_encode('A') | (dna_encode('G') << 2);
    if (c == 'Y')
        return dna_encode('C') | (dna_encode('T') << 2);
    if (c == 'S')
        return dna_encode('G') | (dna_encode('C') << 2);
    if (c == 'W')
        return dna_encode('A') | (dna_encode('T') << 2);
    if (c == 'K')
        return dna_encode('G') | (dna_encode('T') << 2);
    if (c == 'M')
        return dna_encode('A') | (dna_encode('C') << 2);
    if (c == 'B')
        return dna_encode('C') | (dna_encode('G') << 2) | (dna_encode('T') << 4);
    if (c == 'D')
        return dna_encode('A') | (dna_encode('G') << 2) | (dna_encode('T') << 4);
    if (c == 'H')
        return dna_encode('A') | (dna_encode('C') << 2) | (dna_encode('G') << 4);
    if (c == 'N')
        return dna_encode('A') | (dna_encode('C') << 2) | (dna_encode('G') << 4) | (dna_encode('T') << 6);

    return -1;
}

void count_iupac_kmer(uint64_t iupac_kmer, int **dna_cnt)
{
    uint16_t dna_kmer = 0;

    for (int i = 0; i < 8; ++i)
    {
        dna_kmer |= ((iupac_kmer & (1ul << (i * 8))) >> (i * 8));// << 14;
        dna_kmer |= ((iupac_kmer & (1ul << (i * 8 + 1))) >> (i * 8));// << 15;
        if (i != 7)
            dna_kmer <<= 2;
    }

    (*dna_cnt)[dna_kmer] += 1;
}


void process_iupac_kmer(uint64_t iupac_kmer, const int cur_pos, int **dna_cnt)
{
    if (cur_pos == 8)
    {
        count_iupac_kmer(iupac_kmer, dna_cnt);
        return;
    }

    process_iupac_kmer(iupac_kmer, cur_pos + 1, dna_cnt);

    uint8_t cur_letter = iupac_kmer >> (cur_pos * 8);
    const uint64_t nullify_cur = ~(kmer_max << (cur_pos * 8));

    while (cur_letter & extended_bits)
    {
        cur_letter >>= 2;
        process_iupac_kmer((iupac_kmer & nullify_cur) | (cur_letter << (cur_pos * 8)), cur_pos + 1, dna_cnt);
    }
}


int main(void)
{
    int *dna_cnt = calloc(mask_size, sizeof(int));
    uint64_t iupac_kmer = 0;
    char str_ptr;

    char *filename = "/home/alexey/CLionProjects/Solution_c/multiline.example.fasta";
    FILE *fin;

    if ((fin = fopen(filename, "r")) == NULL)
    {
        printf("Error: failed to open file.\n");
        return EXIT_FAILURE;
    }

    char *metadata;
    size_t metadata_len = 0;

    if (getline(&metadata, &metadata_len, fin) == -1)
    {
        printf("Error: failed to read metadata.\n");
        return EXIT_FAILURE;
    }
    free(metadata);

    for (int i = 0; i < 8; ++i)
    {
        str_ptr = fgetc(fin);
        if (isspace(str_ptr))
        {
            --i;
            continue; 
        }

        iupac_kmer <<= 8;
        iupac_kmer |= iupac_encode(str_ptr);
    }
    process_iupac_kmer(iupac_kmer, 0, &dna_cnt);

    const uint64_t left_8_bits_nullifier = (UINT64_MAX >> 8);
    while (!feof(fin))
    {
        str_ptr = fgetc(fin);
        if (str_ptr == EOF)
            break;
        if (isspace(str_ptr))
            continue;

        iupac_kmer &= left_8_bits_nullifier;
        iupac_kmer <<= 8;
        iupac_kmer |= iupac_encode(str_ptr);
        process_iupac_kmer(iupac_kmer, 0, &dna_cnt);
    }

    for (int i = 0; i < mask_size; ++i)
    {
        if (!dna_cnt[i])
            continue;

        char dna_kmer[9];
        dna_kmer_decode(i, dna_kmer);

        printf("%s: %d\n", dna_kmer, dna_cnt[i]);
    }

    free(dna_cnt);
    return EXIT_SUCCESS;
}
