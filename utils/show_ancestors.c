#include "src/kraken_stats.h"
#include <getopt.h>
#include <errno.h>

void print_lvl_name(uint64_t internal_id, Taxonomy const* const tx)
{
	uint64_t name_offset = tx->nodes[internal_id].name_offset;
	uint64_t rank_offset = tx->nodes[internal_id].rank_offset;
	uint64_t ext_id = tx->nodes[internal_id].external_id;
	fprintf(stdout, "%s\t%lu\t%s\n", tx->rank_data+rank_offset, ext_id, tx->name_data+name_offset);
	return;
}

void get_all_parent_names(uint64_t internal_id, Taxonomy const* const tx)
{
	uint64_t name_offset = tx->nodes[internal_id].name_offset;
	if (!strncmp(tx->name_data+name_offset, "root", 4)){
		print_lvl_name(internal_id, tx);
		return;
	} else {
		print_lvl_name(internal_id, tx);
		uint64_t parent_id = tx->nodes[internal_id].parent_id;
		get_all_parent_names(parent_id, tx);
		return;
	}
}

void print_usage(void)
{
    fprintf(stderr, "show_parents -t taxid -d taxo.k2d\n");
    fprintf(stderr, "show_parents --taxid taxid --db taxo.k2d\n");
    fprintf(stderr, "\n");
    return;
}

int main(int argc, char* argv[argc])
{
	static struct option long_opts[] =
    {
        {"db", required_argument, 0, 'd'}
        , {"taxid", required_argument, 0, 't'}
    };
	int opt;
	int l_idx;
	char* db_name = 0;
	char* taxid_str = 0;
	while ((opt = getopt_long(argc, argv, "d:t:", long_opts, &l_idx)) != -1){
        switch (opt) {
            case 'd':
                db_name = strndup(optarg, 1024);
                break;
            case 't':
                taxid_str = strndup(optarg, 1024);
                break;
        }
    }
	if (db_name == 0){
		fprintf(stderr, "taxo.k2d file missing\n");
		print_usage();
		return EXIT_FAILURE;
	}
	if (taxid_str == 0){
		fprintf(stderr, "taxid missing\n");
		print_usage();
		return EXIT_FAILURE;
	}
    Taxonomy* tx = tx_create(db_name);
	uint64_t taxid = strtoul(taxid_str, (void*)0, 10);
	uint64_t int_id = get_internal_id(taxid, tx);
    get_all_parent_names(int_id, tx);
    tx_destroy(tx);
    return EXIT_SUCCESS;
}

