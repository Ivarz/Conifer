#include "src/kraken_stats.h"
#include <getopt.h>
#include <errno.h>

void print_usage(void)
{
    fprintf(stderr, "is_a_parent_of_b taxid1 taxid2 taxo.k2d\n");
    fprintf(stderr, "\n");

    return;
}

int main(int argc, char* argv[argc])
{
    if (argc < 4){
        print_usage();
        return EXIT_FAILURE;
    }
    uint64_t taxid1 = strtoul(argv[1], (void*)0, 10);
    uint64_t taxid2 = strtoul(argv[2], (void*)0, 10);
    Taxonomy* tx = tx_create(argv[3]);
    if (is_a_parent_of_b(taxid1, taxid2, tx)){
        printf("1\n");
    } else {
        printf("0\n");
    }
    tx_destroy(tx);
    return EXIT_SUCCESS;
}
