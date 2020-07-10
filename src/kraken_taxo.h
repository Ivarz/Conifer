#include <stdio.h>                         
#include <string.h>
#include <stdint.h>                  
#include <stdlib.h>                          
#include <stdbool.h>                             
#include "uthash.h"
#define FILE_MAGIC "K2TAXDAT"

typedef struct TaxonomyNode TaxonomyNode;        
struct TaxonomyNode {                            
  uint64_t parent_id;     // Must be lower-numbered node
  uint64_t first_child;   // Must be higher-numbered node
  uint64_t child_count;   // Children of a node are in contiguous block
  uint64_t name_offset;   // Location of name in name data super-string
  uint64_t rank_offset;   // Location of rank in rank data super-string
  uint64_t external_id;   // Taxonomy ID for reporting purposes (usually NCBI)                  
  uint64_t godparent_id;  // Reserved for future use to enable faster traversal
};

typedef struct TaxidMap TaxidMap;
struct TaxidMap {
    uint64_t ext_id;           /* key */
    uint64_t int_id;           /* value */
    UT_hash_handle hh;         /* makes this structure hashable */
};

typedef struct Taxonomy Taxonomy;
struct Taxonomy
{
    int64_t node_count;
    int64_t name_data_len;
    int64_t rank_data_len;
    char* name_data;
    char* rank_data;
    TaxonomyNode* nodes;
    TaxidMap* taxid_map;
};

void txm_add_id(TaxidMap** txm, uint64_t eid, uint64_t iid);

Taxonomy* tx_create(char const* fname);
void tx_destroy(Taxonomy* tx);

void tn_print(TaxonomyNode const* const tn);
uint64_t get_internal_id(uint64_t external_id, Taxonomy const* const tx);
bool is_a_parent_of_b(uint64_t a_ext, uint64_t b_ext, Taxonomy const* const tx);
bool is_a_child_of_b(uint64_t a_ext, uint64_t b_ext, Taxonomy const* const tx);
char const* const tx_taxid_name(uint64_t taxid, Taxonomy const* const tx);
