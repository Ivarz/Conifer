#include "src/kraken_stats.h"
#include <getopt.h>
#include <errno.h>

void print_usage(void)
{
    fprintf(stderr, "taxid_name taxid taxo.k2d\n");
    fprintf(stderr, "\n");

    return;
}

int main(int argc, char* argv[argc])
{
    if (argc < 3){
        print_usage();
        return EXIT_FAILURE;
    }
    uint64_t taxid = strtoul(argv[1], (void*)0, 10);
    Taxonomy* tx = tx_create(argv[2]);
    printf("%s\n",tx_taxid_name(taxid, tx));
    tx_destroy(tx);
    return EXIT_SUCCESS;
}
