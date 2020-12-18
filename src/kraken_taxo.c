#include "src/kraken_taxo.h"
#include <errno.h>
#include "uthash.h"

void tn_print(TaxonomyNode const* const tn)
{
	printf("--------------------\n");
	printf("parent_id %lu\n", tn->parent_id);
	printf("first_child %lu\n", tn->first_child);
	printf("child_count %lu\n", tn->child_count);
	printf("name_offset %lu\n", tn->name_offset);
	printf("rank_offset %lu\n", tn->rank_offset);
	printf("external_id %lu\n", tn->external_id);
	printf("godparent_id %lu\n", tn->godparent_id);
	printf("--------------------\n\n");
	return;
}

uint64_t get_internal_id_binsearch(uint64_t external_id, Taxonomy const* const tx, uint64_t begin, uint64_t end)
{
	if (begin - end < 1024){
		for (uint64_t i=begin; i < end; i++){
			if (tx->nodes[i].external_id == external_id){
				return i;
			}
		}
		return 0;
	} else {
		int middle = (end - begin) / 2;
		if (tx->nodes[middle].external_id > external_id){
			return get_internal_id_binsearch(external_id, tx, begin, middle);
		} else {
			return get_internal_id_binsearch(external_id, tx, middle, end);
		}
	}
	return 0;
}

uint64_t get_internal_id(uint64_t external_id, Taxonomy const* const tx)
{
	TaxidMap* s;
	HASH_FIND_INT(tx->taxid_map, &external_id, s);
	return s == NULL ? 0 : s->int_id;
}

bool is_a_parent_of_b(uint64_t a_ext, uint64_t b_ext, Taxonomy const* const tx)
{
	if (a_ext == 1 && a_ext <= b_ext){
		return true;
	}
	uint64_t b_internal_id = get_internal_id(b_ext, tx);
	TaxonomyNode const* b_node = tx->nodes+b_internal_id;
	while (b_node->parent_id){
		if (b_node->external_id == a_ext){
			return true;
		} else {
			b_node = tx->nodes + b_node->parent_id;
		}
	}
	return false;
}

Taxonomy* tx_create(char const* fname)
{
	FILE* taxo_fh = fopen(fname, "rb");
	if (taxo_fh == NULL){
		printf("Can't open %s\n", fname);
		printf("%s\n", strerror(errno));
		fclose(taxo_fh);
		return 0;
	}
	char magic[strlen(FILE_MAGIC) + 1];
	memset(magic, 0, strlen(FILE_MAGIC) + 1);

	fread(magic, 1, strnlen(FILE_MAGIC,16), taxo_fh);
	if (strcmp(magic, FILE_MAGIC) != 0){
		printf("malformed taxonomy file: %s\n", fname);
		printf("%s\n", magic);
		fclose(taxo_fh);
		return 0;
	}

	Taxonomy* tx = malloc(sizeof(*tx));
	tx->node_count = 0;

	fread(&tx->node_count, 1, sizeof(tx->node_count), taxo_fh);
	fread(&tx->name_data_len, 1, sizeof(tx->name_data_len), taxo_fh);
	fread(&tx->rank_data_len, 1, sizeof(tx->rank_data_len), taxo_fh);

	tx->nodes = malloc(tx->node_count*sizeof(*tx->nodes));
	fread(tx->nodes, sizeof(*tx->nodes), tx->node_count, taxo_fh);

	tx->name_data = malloc(tx->name_data_len*sizeof(*tx->name_data));
	fread(tx->name_data, sizeof(*tx->name_data), tx->name_data_len, taxo_fh);

	tx->rank_data = malloc(tx->rank_data_len*sizeof(*tx->rank_data));
	fread(tx->rank_data, sizeof(*tx->rank_data), tx->rank_data_len, taxo_fh);

	/*printf("%s\n", magic);*/
	/*printf("%ld\n", tx->node_count);*/
	/*printf("%ld\n", tx->name_data_len);*/
	/*printf("%s rank_data_len %ld\n", __func__, tx->rank_data_len);*/

	tx->taxid_map = NULL;
	for (int i=0; i < tx->node_count; i++){
		uint64_t ext_id = tx->nodes[i].external_id;
		/*printf("%s\n", (tx->name_data + (tx->nodes[i].name_offset)));*/
		txm_add_id(&tx->taxid_map, ext_id, i);
		/*tn_print(tx->nodes+i);*/
		/*if (tx->nodes[i].child_count > 0){*/
		/*printf("first_child\n");*/
		/*int first_child = tx->nodes[i].first_child;*/
		/*printf("%s\n", (tx->name_data + (tx->nodes[first_child].name_offset)));*/
		/*tn_print(tx->nodes+first_child);*/
		/*}*/
	}
	fclose(taxo_fh);
	return tx;
}

void tx_destroy(Taxonomy* tx)
{
	free(tx->nodes);
	free(tx->name_data);
	free(tx->rank_data);
	/*free(tx->taxid_map);*/
	free(tx);
	return;
}


char* tx_taxid_name(uint64_t taxid, Taxonomy const* const tx)
{
	uint64_t int_id = get_internal_id(taxid, tx);
	uint64_t name_offset = tx->nodes[int_id].name_offset;
	return tx->name_data+name_offset;
}

void txm_add_id(TaxidMap** txm, uint64_t eid, uint64_t iid)
{
	TaxidMap *s;
	HASH_FIND_INT(*txm, &eid, s);  /* id already in the hash? */
	if (s==NULL) {
		s = malloc(sizeof(*s));
		s->ext_id = eid;
		s->int_id = iid;
		HASH_ADD_INT( *txm, ext_id, s);  /* id: name of key field */
	}
	return;
}
