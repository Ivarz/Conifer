#ifndef KRAKEN_STATS_H
#define KRAKEN_STATS_H
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <limits.h>
#include "src/kraken_taxo.h"

#define LINE_SIZE 4096
#define KMER_COUNT_SIZE 1024
#define TAXID_KMER_DATA_SIZE 1024

inline
void check_ulong_overflow(uint64_t value)
{
    if (value == ULLONG_MAX && errno == ERANGE){
        fprintf(stderr, "CRITICAL_ERROR: value bigger than %lu\n", ULONG_MAX);
        exit(1);
    }
    return;
}


typedef struct KmerCounts KmerCounts;
struct KmerCounts
{
    size_t size;
    size_t capacity;
    uint64_t* taxids;
    uint64_t* counts;
};

KmerCounts* kmc_create(void);
KmerCounts* kmc_reset(KmerCounts* kc);
KmerCounts* kmc_fill(KmerCounts* kc, char* str);
void kmc_destroy(KmerCounts* kc);

typedef struct KrakenRec KrakenRec;
struct KrakenRec
{
    bool paired;
    bool classified;
    char* read_name;
    uint64_t taxid;
    uint64_t read1_len;
    uint64_t read2_len;
    KmerCounts* read1_kmers;
    KmerCounts* read2_kmers;
};

KrakenRec* kraken_create(bool paired);
KrakenRec* kraken_fill(KrakenRec* krp, char* kraken_line);
void kraken_print(KrakenRec* krp);
KrakenRec* kraken_reset(KrakenRec* krp);
void kraken_destroy(KrakenRec* krp);

typedef struct FloatHist FloatHist;
struct FloatHist
{
    size_t size;
    size_t capacity;
    int32_t* fraction_counts;
};

FloatHist* fh_create();
void fh_destroy(FloatHist* fh);
void fh_add(FloatHist* fh, float frac);
void fh_print(FloatHist* fh);

typedef struct TaxIdData TaxIdData;
struct TaxIdData
{
    int32_t taxid_size;
    int32_t taxid_capacity;
    uint64_t* taxids;
    FloatHist** data;
};

TaxIdData* txd_create();
void txd_destroy(TaxIdData* txd);
int32_t txd_get_taxa_index(TaxIdData* txd, uint64_t taxid);
TaxIdData* txd_extend(TaxIdData* txd, int32_t new_capacity);
int32_t txd_add_new_taxa(TaxIdData* txd, uint64_t const taxid);
void txd_add_data(TaxIdData* txd, uint64_t taxid, float frac);
void txd_print(TaxIdData* txd);

typedef struct Quartiles Quartiles;
struct Quartiles
{
    float q1;
    float q2;
    float q3;
};

int64_t fh_sum(FloatHist* fh);
Quartiles get_quartiles(FloatHist* fh);
float kmer_fraction(KmerCounts const* const kmcs, uint64_t const taxid);

float get_avg_kmer_fraction(KrakenRec* krp);
KrakenRec* kraken_adjust_taxonomy(KrakenRec* krp, Taxonomy const* const tx);
KrakenRec* kraken_adjust_taxonomy_nonconflicting(KrakenRec* krp, Taxonomy const* const tx);

#endif
