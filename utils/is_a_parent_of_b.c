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

void print_result(uint64_t t1, uint64_t t2, bool is_related_flag, bool labels_flag, Taxonomy const* tx)
{
    if (labels_flag){
        printf("%lu\t%lu\t", t1, t2);
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
        , {"taxid1", required_argument, 0, '1'}
        , {"taxid2", required_argument, 0, '2'}
        , {"any", no_argument, 0, 'a'}
        , {"labels", no_argument, 0, 'b'}
    };
    int opt;
    if (argc < 4){
        print_usage();
        return EXIT_FAILURE;
    }
    char* db_name = 0;
    char* file_name = 0;
    uint64_t taxid1 = 0;
    uint64_t taxid2 = 0;
    int l_idx = 0;
    bool is_related_flag = false;
    bool labels_flag = false;
    while ((opt = getopt_long(argc, argv, "d:l:1:2rb", long_opts, &l_idx)) != -1){
        switch (opt) {
            case 'l':
                file_name = strndup(optarg, 1024);
                break;
            case 'd':
                db_name = strndup(optarg, 1024);
                break;
            case '1':
                taxid1 = strtoul(optarg, (void*)0, 10);
                break;
            case '2':
                taxid2 = strtoul(optarg, (void*)0, 10);
                break;
            case 'a':
                is_related_flag = true;
                break;
            case 'b':
                labels_flag = true;
                break;
        }
    }
    Taxonomy* tx = tx_create(db_name);
    if (!file_name){
        print_result(taxid1, taxid2, is_related_flag, labels_flag, tx);
    } else {
        char line[LINE_SIZE];
        FILE* fh = fopen(file_name, "r");
        while(fgets(line, sizeof(line), fh)){
            uint64_t taxid1 = strtoul(strtok(line,"\t"), (void*)0, 10);
            uint64_t taxid2 = strtoul(strtok(0,"\n"), (void*)0, 10);
            print_result(taxid1, taxid2, is_related_flag, labels_flag, tx);
        }
        fclose(fh);
        free(file_name);
    }

    tx_destroy(tx);
    free(db_name);
    return EXIT_SUCCESS;
}
