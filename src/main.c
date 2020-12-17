#include "src/kraken_stats.h"
#include "src/error_type.h"
#include <getopt.h>
#include <errno.h>
#include <assert.h>
#include <zlib.h>

#define MAJOR 1
#define MINOR 0
#define PATCH 2

#define LINE_SIZE 4096
#define SUMMARY 1
#define RTL (1 << 1)
#define BOTH_SCORES (1 << 2)
#define ALL_RECORDS (1 << 3)
#define FILTER (1 << 4)
#define SHOW_VERSION (1 << 5)

void print_usage(void)
{
    fprintf(stderr, "Conifer %d.%d.%d\n", MAJOR, MINOR, PATCH);
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "conifer [OPTIONS] -i <KRAKEN_FILE> -d <TAXO_K2D>\n");
    fprintf(stderr, "\t-i,--input\t\tinput file\n");
    fprintf(stderr, "\t-d,--db\t\t\tkraken2 taxo.k2d file\n");
    fprintf(stderr, "\t-a,--all\t\toutput all reads (including unclassified)\n");
    fprintf(stderr, "\t-s,--summary\t\toutput summary statistics for each taxonomy\n");
    fprintf(stderr, "\t-f,--filter\t\tfilter kraken file by confidence score\n");
    fprintf(stderr, "\t-r,--rtl\t\treport root-to-leaf score instead of confidence score\n");
    fprintf(stderr, "\t-b,--both_scores\treport confidence and root-to-leaf score\n");
    fprintf(stderr, "\t-v,--version\t\tshow version\n");
    fprintf(stderr, "\n");

    return;
}

void print_output(char const* const line, size_t const kmfn, KmerFractions kmfs[kmfn])
{
    printf("%s", line);
    if (kmfs[0].paired){
        for (size_t i = 0; i < kmfn; i++){
            printf("\t%.4f\t%.4f\t%.4f"
                    , kmfs[i].read1_kmer_frac
                    , kmfs[i].read2_kmer_frac
                    , kmfs[i].avg_kmer_frac
                    );
        }
        printf("\n");
    } else {
        for (size_t i = 0; i < kmfn; i++){
            printf("\t%.4f"
                    , kmfs[i].avg_kmer_frac
                    );
        }
        printf("\n");
    }
    return;
}

ErrorType gather_and_print_summary(gzFile fh, Taxonomy const* const tx, int flags)
{
    /*char line[LINE_SIZE] = {0};*/
    /*char line_to_parse[LINE_SIZE] = {0};*/
	String* line = string_create();
	String* line_cpy = string_create();

    int counter = 0;
    int const kinds_of_calculations = (flags & BOTH_SCORES) ? 2 : 1;

    int indices[kinds_of_calculations];

    if (flags & BOTH_SCORES){
        indices[0] = 0;
        indices[1] = 1;
    } else if (flags & RTL){
        indices[0] = 1;
    } else {
        indices[0] = 0;
    }

    KrakenRec* krp = kraken_create(true);
    TaxIdData* txds[2] = {txd_create(), txd_create()};

	KrakenRec* (*tax_adj_fp[2])(KrakenRec*, Taxonomy const* const) = {kraken_adjust_taxonomy, kraken_adjust_taxonomy_rtl};
	ErrorType parsing_status = Success;
    while (parse_line(fh, line) && parsing_status == Success){
		for (int i = 0; i < kinds_of_calculations; i++){
			string_copy(line_cpy, line);
			int j = indices[i];
			parsing_status = kraken_fill(krp, line_cpy);
			if (parsing_status != Success){
				break;
			}
			if (krp->taxid > 0){
				krp = tax_adj_fp[j](krp, tx);
				KmerFractions kmf = kmf_calculate(krp);
				txd_add_data(txds[i], krp->taxid, kmf.avg_kmer_frac);
			}
			krp = kraken_reset(krp);
		}
		if (!(counter % 1000000) && counter){
			fprintf(stderr, "\r%d lines processed...", counter);
			fflush(stderr);
		}
		counter++;
		string_reset(line);
    }
    if (flags & BOTH_SCORES){
        printf("taxon_name\ttaxid\treads\tP25_conf\tP50_conf\tP75_conf\tP25_rtl\tP50_rtl\tP75_rtl\n");
        for (int i=0; i < txds[0]->taxid_size; i++){
            Quartiles qs1 = get_quartiles(txds[0]->data[i]);
            Quartiles qs2 = get_quartiles(txds[1]->data[i]);
            assert(txds[0]->taxids[i] == txds[1]->taxids[i]);
            assert(fh_sum(txds[0]->data[i]) == fh_sum(txds[1]->data[i]));
            printf("%s\t%lu\t%ld\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n"
                    , tx_taxid_name(txds[0]->taxids[i], tx)
                    , txds[0]->taxids[i]
                    , fh_sum(txds[0]->data[i])
                    , qs1.q1/1000.0f, qs1.q2/1000.0f, qs1.q3/1000.0f
                    , qs2.q1/1000.0f, qs2.q2/1000.0f, qs2.q3/1000.0f
                    );
        }
    } else {
        printf("taxon_name\ttaxid\treads\tP25\tP50\tP75\n");
        for (int i=0; i < txds[0]->taxid_size; i++){
            Quartiles qs = get_quartiles(txds[0]->data[i]);
            printf("%s\t%lu\t%ld\t%.4f\t%.4f\t%.4f\n"
                    , tx_taxid_name(txds[0]->taxids[i], tx)
                    , txds[0]->taxids[i]
                    , fh_sum(txds[0]->data[i])
                    , qs.q1/1000.0f, qs.q2/1000.0f, qs.q3/1000.0f);
        }
    }
    txd_destroy(txds[0]);
    txd_destroy(txds[1]);
	kraken_destroy(krp);
	string_destroy(line_cpy);
	string_destroy(line);

    return parsing_status;
}

ErrorType print_scores_by_record(gzFile fh, Taxonomy const* const tx, int flags, float filter_threshold)
{
    /*char line[LINE_SIZE] = {0};*/
    /*char line_to_parse[LINE_SIZE] = {0};*/
	String* line = string_create();
	String* line_cpy = string_create();
    int counter = 1;

    int const kinds_of_calculations = (flags & BOTH_SCORES) ? 2 : 1;
    int indices[kinds_of_calculations];

    if (flags & BOTH_SCORES){
        indices[0] = 0;
        indices[1] = 1;
    } else if (flags & RTL){
        indices[0] = 1;
    } else {
        indices[0] = 0;
    }

    KrakenRec* krp = kraken_create(true);
    KmerFractions kmfs[2];

    KrakenRec* (*tax_adj_fp[2])(KrakenRec*, Taxonomy const* const) = {kraken_adjust_taxonomy, kraken_adjust_taxonomy_rtl};
	ErrorType parsing_status = Success;
    while (parse_line(fh, line) && parsing_status == Success){
        for (int i = 0; i < kinds_of_calculations; i++){
			string_copy(line_cpy, line);
            int j = indices[i];
			parsing_status = kraken_fill(krp, line_cpy);
			if (parsing_status != Success){
				fprintf(stderr, "Malformed input at line %d\n", counter);
				break;
			}
            if (krp->taxid > 0){
                krp = tax_adj_fp[j](krp, tx);
                kmfs[i] = kmf_calculate(krp);
            }
        }
        if (krp->taxid > 0){
            if ((flags & FILTER) && !(flags & BOTH_SCORES)){
                if (kmfs[0].avg_kmer_frac >= filter_threshold){
                    print_output(line->str, kinds_of_calculations, kmfs);
                }
            } else {
                print_output(line->str, kinds_of_calculations, kmfs);
            }
        } else if(flags & ALL_RECORDS){
            kmfs[0].paired = krp->paired;
            kmfs[0].read1_kmer_frac = 0.0f;
            kmfs[0].read2_kmer_frac = 0.0f;
            kmfs[0].avg_kmer_frac = 0.0f;
            kmfs[1] = kmfs[0];
            print_output(line->str, kinds_of_calculations, kmfs);
        }
        if (!(counter % 1000000) && counter){
            fprintf(stderr, "\r%d lines processed...", counter);
			fflush(stderr);
        }
        counter++;
        krp = kraken_reset(krp);
		string_reset(line);
    }

	string_destroy(line_cpy);
	string_destroy(line);
	kraken_destroy(krp);
    return parsing_status;
}

int main(int argc, char* argv[argc])
{
    static struct option long_opts[] =
    {
        {"input", required_argument, 0, 'i'}
        , {"db", required_argument, 0, 'd'}
        , {"all_reads", no_argument, 0, 'a'}
        , {"summary", no_argument, 0, 's'}
        , {"rtl", no_argument, 0, 'r'}
        , {"filter", required_argument, 0, 'f'}
        , {"both_scores", no_argument, 0, 'b'}
        , {"version", no_argument, 0, 'v'}

    };
    if (argc < 2){
        print_usage();
        return EXIT_FAILURE;
    }
    int opt;
    char* file_name = 0;
    char* db_name = 0;
    int flags = 0;
    int l_idx = 0;
    float filter_threshold = -1.0f;
    char* filter_threshold_str = 0;
    while ((opt = getopt_long(argc, argv, "i:d:rbaspv", long_opts, &l_idx)) != -1){
        switch (opt) {
            case 'i':
                file_name = strndup(optarg, 1024);
                break;
            case 'd':
                db_name = strndup(optarg, 1024);
                break;
            case 'a':
                flags |= ALL_RECORDS;
                break;
            case 's':
                flags |= SUMMARY;
                break;
            case 'r':
                flags |= RTL;
                break;
            case 'b':
                flags |= BOTH_SCORES;
                break;
            case 'v':
                flags |= SHOW_VERSION;
                break;
            case 'f':
                filter_threshold_str = strndup(optarg, 1024);
                flags |= FILTER;
                break;
        }
    }
	if (flags & SHOW_VERSION){
		fprintf(stderr, "Conifer %d.%d.%d\n", MAJOR, MINOR, PATCH);
        return EXIT_SUCCESS;
	}
    if ((flags & ALL_RECORDS) && (flags & SUMMARY)){
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
    if (flags & FILTER){
        filter_threshold = strtod(filter_threshold_str, 0);
        free(filter_threshold_str);
    }

    gzFile fh = gzopen(file_name, "rb");

    if (fh == NULL){
        fprintf(stderr, "Can't open %s\n", file_name);
        fprintf(stderr, "%s\n", strerror(errno));
        return EXIT_FAILURE;
    }

    Taxonomy* tx = tx_create(db_name);
    if (tx == NULL){
        return EXIT_FAILURE;
    }

	int exit_status = EXIT_FAILURE;
    if (flags & SUMMARY){
        if (gather_and_print_summary(fh, tx, flags)){
			exit_status = EXIT_SUCCESS;
		}
    } else {
        if(print_scores_by_record(fh, tx, flags, filter_threshold)){
			exit_status = EXIT_SUCCESS;
		}
    }
    tx_destroy(tx);
    gzclose(fh);
    free(file_name);
    free(db_name);
    return exit_status;
}
