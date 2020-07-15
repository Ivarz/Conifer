#include "src/kraken_stats.h"
#include <getopt.h>
#include <errno.h>

void print_usage(void)
{
    fprintf(stderr, "conifer [OPTIONS] -i <KRAKEN_FILE> -d <TAXO_K2D>\n");
    fprintf(stderr, "\t-i,--input\t\tinput file\n");
    fprintf(stderr, "\t-d,--db\t\tkraken2 taxo.k2d file\n");
    fprintf(stderr, "\t-a,--all\t\toutput all reads (including unclassified)\n");
    fprintf(stderr, "\t-s,--summary\t\toutput summary statistics for each taxonomy\n");
    fprintf(stderr, "\n");

    return;
}

int main(int argc, char* argv[argc])
{
    static struct option long_opts[] =
    {
        {"input", required_argument, 0, 'i'}
        , {"db", required_argument, 0, 'd'}
        , {"all", no_argument, 0, 'a'}
        , {"summary", no_argument, 0, 's'}
        , {"non_conflicting", no_argument, 0, 'n'}

    };
    if (argc < 3){
        print_usage();
        return EXIT_FAILURE;
    }
    int opt;
    char* file_name = 0;
    char* db_name = 0;
    bool summary = false;
    bool all_reads = false;
    bool non_conflicting = false;
    int l_idx = 0;
    while ((opt = getopt_long(argc, argv, "i:d:nbasp", long_opts, &l_idx)) != -1){
        switch (opt) {
            case 'i':
                file_name = strndup(optarg, 1024);
                break;
            case 'd':
                db_name = strndup(optarg, 1024);
                break;
            case 'a':
                all_reads = true;
                break;
            case 's':
                summary = true;
                break;
            case 'n':
                non_conflicting = true;
                break;
        }
    }
    if (all_reads && summary){
        fprintf(stderr, "Options -a or -s are mutually exclusive\n");
        return EXIT_FAILURE;
    }
    if (file_name == NULL){
        fprintf(stderr, "Provide input file name \n");
        return EXIT_FAILURE;
    }
    if (db_name == NULL){
        fprintf(stderr, "Provide kraken2 taxo.k2d filename\n");
        return EXIT_FAILURE;
    }

    FILE* fh = fopen(file_name, "r");

    if (fh == NULL){
        fprintf(stderr, "Can't open %s\n", file_name);
        fprintf(stderr, "%s\n", strerror(errno));
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

    KrakenRec* (*tax_adj_fp)(KrakenRec*, Taxonomy const* const) = non_conflicting ? kraken_adjust_taxonomy_nonconflicting : kraken_adjust_taxonomy;
    while(fgets(line, sizeof(line), fh)){
        memcpy(line_to_parse, line, sizeof(*line)*LINE_SIZE);
        krp = kraken_fill(krp, line_to_parse);
        if (krp->taxid > 0){
            krp = tax_adj_fp(krp, tx);
            float avg_kmer_frac = -1.0f;
            float kmer_frac1 = -1.0f;
            float kmer_frac2 = -1.0f;
            if (krp->paired){
                kmer_frac1 = krp->read1_kmers->size ? kmer_fraction(krp->read1_kmers, krp->taxid) : -1.0f;
                kmer_frac2 = krp->read2_kmers->size ? kmer_fraction(krp->read2_kmers, krp->taxid) : -1.0f;
                avg_kmer_frac = kmer_frac1 == -1.0f && kmer_frac2 == -1.0f ? 0.0f :
                    (kmer_frac1 >= 0.0f && kmer_frac2 < 0.0f) ? kmer_frac1 :
                    (kmer_frac1 < 0.0f && kmer_frac2 >= 0.0f) ? kmer_frac2 :
                    ((kmer_frac1 + kmer_frac2) / 2.0f);
            } else {
                kmer_frac1 = krp->read1_kmers->size ? kmer_fraction(krp->read1_kmers, krp->taxid) : -1.0f;
                avg_kmer_frac = kmer_frac1 == -1.0f ? 0.0f : kmer_frac1;
            }
            if (summary){
                txd_add_data(txd, krp->taxid, avg_kmer_frac);
            } else {
                line[strnlen(line, LINE_SIZE) -1] = '\0';
                if (krp->paired){
                    printf("%s\t%.4f\t%.4f\t%.4f\n", line, kmer_frac1, kmer_frac2, avg_kmer_frac);
                } else {
                    printf("%s\t%.4f\n", line, avg_kmer_frac);
                }
            }
        } else {
            if(all_reads){
                line[strnlen(line, LINE_SIZE) -1] = '\0';
                if (krp->paired){
                    printf("%s\t%.4f\t%.4f\t%.4f\n", line, 0.0f, 0.0f, 0.0f);
                } else {
                    printf("%s\t%.4f\n", line, 0.0f);
                }

            }
        }
        krp = kraken_reset(krp);
        if (!(counter % 1000000) && counter){
            fprintf(stderr, "%d lines processed...\n", counter);
        }
        counter++;
    }

    if (summary)
    printf("taxon_name\ttaxid\treads\tP25\tP50\tP75\n");
    for (int i=0; i < txd->taxid_size; i++){
        Quartiles qs = get_quartiles(txd->data[i]);
        printf("%s\t%lu\t%ld\t%.4f\t%.4f\t%.4f\n"
                , tx_taxid_name(txd->taxids[i], tx)
                , txd->taxids[i]
                , fh_sum(txd->data[i])
                , qs.q1/1000.0f, qs.q2/1000.0f, qs.q3/1000.0f);
    }

    kraken_destroy(krp);
    txd_destroy(txd);
    tx_destroy(tx);
    fclose(fh);
    free(file_name);
    free(db_name);
    return EXIT_SUCCESS;
}
