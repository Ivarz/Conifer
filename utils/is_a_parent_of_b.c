#include "src/kraken_stats.h"
#include <getopt.h>
#include <errno.h>

#define LINE_SIZE 4096

void print_usage(void)
{
    fprintf(stderr, "is_a_parent_of_b -1 taxid1 -2 taxid2 --db taxo.k2d\n");
    fprintf(stderr, "is_a_parent_of_b --list taxid_pair_list --db taxo.k2d\n");
    fprintf(stderr, "is_a_parent_of_b --any --labels --list taxid_pair_list --db taxo.k2d\n");
    fprintf(stderr, "\n");

    return;
}


bool is_related(uint64_t t1, uint64_t t2, Taxonomy const* tx)
{
    bool t1_parent_t2 = is_a_parent_of_b(t1, t2, tx);
    bool t2_parent_t1 = is_a_parent_of_b(t2, t1, tx);
    return t1_parent_t2 || t2_parent_t1;
}

void print_result(uint64_t t1, uint64_t t2, bool is_related_flag, bool labels_flag, bool names_flag, Taxonomy const* tx)
{
    if (labels_flag){
        printf("%lu\t%lu\t", t1, t2);
    }
    if (names_flag){
        printf("%s\t%s\t", tx_taxid_name(t1, tx), tx_taxid_name(t2, tx));
    }
    if (is_related_flag){
        if (is_related(t1, t2, tx)){
            printf("1\n");
        } else {
            printf("0\n");
        }
    } else {
        bool t1_parent_t2 = is_a_parent_of_b(t1, t2, tx);
        if (t1_parent_t2){
            printf("1\n");
        } else {
            printf("0\n");
        }
    }
    return;
}

int main(int argc, char* argv[argc])
{
    static struct option long_opts[] =
    {
        {"db", required_argument, 0, 'd'}
        , {"list", required_argument, 0, 'l'}
		, {"csv", required_argument, 0, 'c'}
        , {"taxid1", required_argument, 0, '1'}
        , {"taxid2", required_argument, 0, '2'}
        , {"any", no_argument, 0, 'a'}
        , {"labels", no_argument, 0, 'b'}
        , {"names", no_argument, 0, 'n'}
    };
    int opt;
    char* db_name = 0;
    char* file_name = 0;
    char* csv_name = 0;
    char* taxid1_str = 0;
    char* taxid2_str = 0;
    uint64_t taxid1 = 0;
    uint64_t taxid2 = 0;
    int l_idx = 0;
    bool is_related_flag = false;
    bool labels_flag = false;
    bool names_flag = false;
    while ((opt = getopt_long(argc, argv, "d:l:c:1:2:rbn", long_opts, &l_idx)) != -1){
        switch (opt) {
            case 'l':
                file_name = strndup(optarg, 1024);
                break;
            case 'c':
                csv_name = strndup(optarg, 1024);
                break;
            case 'd':
                db_name = strndup(optarg, 1024);
                break;
            case '1':
                taxid1_str = strndup(optarg, 1024);
                break;
            case '2':
                taxid2_str = strndup(optarg, 1024);
                break;
            case 'a':
                is_related_flag = true;
                break;
            case 'b':
                labels_flag = true;
                break;
            case 'n':
                names_flag = true;
                break;
        }
    }
    if (argc < 4 || !db_name){
        print_usage();
        return EXIT_FAILURE;
    }
    Taxonomy* tx = tx_create(db_name);
    if (file_name) {
        char line[LINE_SIZE];
        FILE* fh = fopen(file_name, "r");
        while(fgets(line, sizeof(line), fh)){
            uint64_t taxid1 = strtoul(strtok(line,"\t"), (void*)0, 10);
            uint64_t taxid2 = strtoul(strtok(0,"\n"), (void*)0, 10);
			print_result(taxid1, taxid2, is_related_flag, labels_flag, names_flag, tx);
        }
        fclose(fh);
        free(file_name);
    } else if (csv_name){
        char line[LINE_SIZE];
        char* line_cpy;
        FILE* fh = fopen(csv_name, "r");
        while(fgets(line, sizeof(line), fh)){
			line_cpy = strndup(line, LINE_SIZE);
			char* taxid1_str = strtok(line,"\t");
			char* taxid2_str = strtok(0,",\n");
            uint64_t taxid1 = strtoul(taxid1_str, (void*)0, 10);
			bool is_hit = false;
			while (taxid2_str){
				/*printf("%s\n", taxid2_str);*/
				uint64_t taxid2 = strtoul(taxid2_str, (void*)0, 10);
				is_hit = is_hit ? is_hit : is_related(taxid1, taxid2, tx);
				taxid2_str = strtok(0,",\n");
			}
			line_cpy[strnlen(line_cpy, LINE_SIZE) -1] = '\0';
			printf("%s\t%d\n", line_cpy, (int) is_hit);
			free(line_cpy);
        }
        fclose(fh);
        free(csv_name);
	} else {
        /*printf("%s\t%s\n", taxid1_str, taxid2_str);*/
        taxid1 = strtoul(taxid1_str, (void*)0, 10);
        taxid2 = strtoul(taxid2_str, (void*)0, 10);
        print_result(taxid1, taxid2, is_related_flag, labels_flag, names_flag, tx);
        free(taxid1_str);
        free(taxid2_str);
	}

    tx_destroy(tx);
    free(db_name);
    return EXIT_SUCCESS;
}
