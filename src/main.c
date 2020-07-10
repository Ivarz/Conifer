#include "src/kraken_stats.h"
#include <getopt.h>
#include <errno.h>

int main(int argc, char* argv[argc])
{
    int opt;
    char* file_name = 0;
    char* db_name = 0;
    while ((opt = getopt(argc, argv, "i:d:")) != -1){
        switch (opt) {
            case 'i':
                file_name = strndup(optarg, 1024);
                break;
            case 'd':
                db_name = strndup(optarg, 1024);
                break;
        }
    }
    if (file_name == NULL){
        printf("Provide input file name \n");
        return EXIT_FAILURE;
    }
    if (db_name == NULL){
        printf("Provide kraken2 taxo.k2d filename\n");
        return EXIT_FAILURE;
    }

    FILE* fh = fopen(file_name, "r");

    if (fh == NULL){
        printf("Can't open %s\n", file_name);
        printf("%s\n", strerror(errno));
        return EXIT_FAILURE;
    }
    
    Taxonomy* tx = tx_create(db_name);
    if (tx == NULL){
        return EXIT_FAILURE;
    }
    char line[LINE_SIZE] = {0};
    char line_to_parse[LINE_SIZE] = {0};
    int counter = 0;
    KrakenRec* krp = kraken_create(true);
    TaxIdData* txd = txd_create();
    bool summary = false;
    bool all_summary = false;
    while(fgets(line, sizeof(line), fh)){
        memcpy(line_to_parse, line, sizeof(*line)*LINE_SIZE);
        krp = kraken_fill(krp, line_to_parse);
        if (krp->taxid > 0){
            krp = kraken_adjust_taxonomy(krp, tx);
            float avg_kmer_frac = -1.0f;
            if (krp->paired){
                float kmer_frac1 = krp->read1_kmers->size ? kmer_fraction(krp->read1_kmers, krp->taxid) : -1.0f;
                float kmer_frac2 = krp->read2_kmers->size ? kmer_fraction(krp->read2_kmers, krp->taxid) : -1.0f;
                avg_kmer_frac = kmer_frac1 == -1.0f && kmer_frac2 == -1.0f ? 0.0f :
                    (kmer_frac1 >= 0.0f && kmer_frac2 < 0.0f) ? kmer_frac1 :
                    (kmer_frac1 < 0.0f && kmer_frac2 >= 0.0f) ? kmer_frac2 :
                    ((kmer_frac1 + kmer_frac2) / 2.0f);
            } else {
                float kmer_frac1 = krp->read1_kmers->size ? kmer_fraction(krp->read1_kmers, krp->taxid) : -1.0f;
                avg_kmer_frac = kmer_frac1 == -1.0f ? 0.0f : kmer_frac1;
            }
            if (summary){
                txd_add_data(txd, krp->taxid, avg_kmer_frac);
            } else {
                line[strnlen(line, LINE_SIZE) -1] = '\0';
                printf("%s\t%.4f\n", line, avg_kmer_frac);
            }
        } else {
            if(all_summary){
                line[strnlen(line, LINE_SIZE) -1] = '\0';
                printf("%s\t%.4f\n", line, 0.0f);
            }
        }
        krp = kraken_reset(krp);
        if (!(counter % 1000000) && counter){
            fprintf(stderr, "%d lines processed...\n", counter);
        }
        counter++;
    }

    if (summary)
    printf("taxon_name\ttaxid\tQ1\tQ2\tQ3\n");
    for (int i=0; i < txd->taxid_size; i++){
        Quartiles qs = get_quartiles(txd->data[i]);
        printf("%s\t%lu\t%.4f\t%.4f\t%.4f\n", tx_taxid_name(txd->taxids[i], tx), txd->taxids[i], qs.q1/1000.0f, qs.q2/1000.0f, qs.q3/1000.0f);
    }

    kraken_destroy(krp);
    txd_destroy(txd);
    tx_destroy(tx);
    fclose(fh);
    free(file_name);
    free(db_name);
    return EXIT_SUCCESS;
}
