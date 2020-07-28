#include "src/kraken_stats.h"

void check_ulong_overflow(uint64_t);

KmerCounts* kmc_create(void)
{
    KmerCounts* kc = malloc(sizeof(*kc));
    if (kc == 0){
        printf("Failed to malloc kmer counts\n");
        return 0;
    }

    kc->size = 0;
    kc->capacity = KMER_COUNT_SIZE;
    kc->taxids = malloc(sizeof(*kc->taxids) * kc->capacity);
    kc->counts = malloc(sizeof(*kc->counts) * kc->capacity);
    if (kc->taxids == 0 || kc->counts == 0){
        free(kc);
        printf("Failed to malloc kmer counts\n");
        return 0;
    }
    return kc;
}

KmerCounts* kmc_reset(KmerCounts* kc)
{
    kc->size = 0;
    memset(kc->taxids, 0, sizeof(*kc->taxids)*kc->capacity);
    memset(kc->counts, 0, sizeof(*kc->counts)*kc->capacity);
    return kc;
}

KmerCounts* kmc_extend(KmerCounts* kc, size_t new_size)
{
    kc->taxids = realloc(kc->taxids, sizeof(*kc->taxids) * new_size);
    kc->counts = realloc(kc->counts, sizeof(*kc->counts) * new_size);
    if (kc->taxids == 0 || kc->counts == 0){
        printf("Failed to malloc kmer counts\n");
        /*free(kc);*/
        return 0;
    }
    return kc;
}


KmerCounts* kmc_fill(KmerCounts* kc, char* str)
{
    if (!strcmp(str, "") || !strcmp(str,":")){
        kc = kmc_reset(kc);
        return kc;
    }
    char* taxid_str = strtok(str,":");
    char* count_str = strtok(0," ");
    size_t i = 0;
    kc->size = 0;

    while(taxid_str != 0 && count_str != 0){
        if (i >= kc->capacity){
            kc->capacity += KMER_COUNT_SIZE;
            kc = kmc_extend(kc, kc->capacity);
            if (kc == 0){
                fprintf(stderr, "Failed to extend kmc at line %d", __LINE__);
            }
        }
        if (taxid_str[0] != 'A'){
            kc->taxids[i] = strtoul(taxid_str, (void*)0, 10);
            check_ulong_overflow(kc->taxids[i]);
            kc->counts[i] = strtoul(count_str, (void*)0, 10);
            check_ulong_overflow(kc->counts[i]);
            i++;
        }
        taxid_str = strtok(0,":");
        count_str = strtok(0," ");
    }
    kc->size = i;
    return kc;
}

/*void kmc_copy(KmerCounts const* const src, KmerCounts* const dst)*/
/*{*/
    /*return;*/
/*}*/
void kmc_destroy(KmerCounts* kc)
{
    free(kc->counts);
    free(kc->taxids);
    free(kc);
    return;
}

KrakenRec* kraken_create(bool paired)
{
    KrakenRec* krp = malloc(sizeof(*krp));
    krp->paired = paired;
    krp->classified = false;
    krp->read_name = "";
    krp->taxid = 0;
    krp->read1_len = 0;
    krp->read2_len = 0;
    krp->read1_kmers = kmc_create();
    krp->read2_kmers = kmc_create();
    return krp;
}

KrakenRec* kraken_reset(KrakenRec* krp)
{
    krp->classified = false;
    krp->read_name = "";
    krp->taxid = 0;
    krp->read1_len = 0;
    krp->read2_len = 0;
    krp->read1_kmers = kmc_reset(krp->read1_kmers);
    krp->read2_kmers = kmc_reset(krp->read2_kmers);
    return krp;
}

bool paired_end_data(char* kraken_line)
{
    bool paired = false;
    for (size_t i = 0; i < strnlen(kraken_line, 1024); i++){
        if (kraken_line[i] == '|'){
            paired = true;
            return paired;
        }
    }
    return paired;
}

KrakenRec* kraken_fill(KrakenRec* krp, char* kraken_line)
{
    if (!strcmp(kraken_line,"")){
        return kraken_reset(krp);
    }
    if (kraken_line[0] == 'U'){
        return kraken_reset(krp);
    }
    krp->paired = paired_end_data(kraken_line);
    if (krp->paired){
        char* class_field = strtok(kraken_line, "\t");
        char* read_name = strtok(0, "\t");
        char* tax_field = strtok(0, "\t");
        char* read1_len = strtok(0, "|");
        char* read2_len = strtok(0, "\t");
        char* kmer_pairs_str = strtok(0, "\t");
        char* kmer_pairs_str1 = strtok(kmer_pairs_str, "|");
        //if missing kmer_pairs_str1 adjust offest for second read kmers
        int offset =
            !strcmp(kmer_pairs_str1,"") ? 3 :
            !strcmp(kmer_pairs_str1,":") ? 1 : 3;
        char* kmer_pairs_str2 = strtok(0, "\n") + offset;
        /*printf("%s\n", class_field);*/
        /*printf("%s\n", read_name);*/
        /*printf("%s\n", tax_field);*/
        /*printf("%s\n", read1_len);*/
        /*printf("%s\n", read2_len);*/
        /*printf("%s\n", kmer_pairs_str1);*/
        /*printf("%s\n", kmer_pairs_str2);*/

        krp->classified = class_field[0] == 'C' ? true : false;
        krp->read_name = read_name;
        krp->taxid = strtoul(tax_field, (void*)0, 10);
        check_ulong_overflow(krp->taxid);

        krp->read1_len = strtoul(read1_len, (void*)0, 10);
        check_ulong_overflow(krp->read1_len);
        krp->read2_len = strtoul(read2_len, (void*)0, 10);
        check_ulong_overflow(krp->read2_len);


        krp->read1_kmers = kmc_fill(krp->read1_kmers, kmer_pairs_str1);
        krp->read2_kmers = kmc_fill(krp->read2_kmers, kmer_pairs_str2);
        return krp;
    } else {
        char* class_field = strtok(kraken_line, "\t");
        char* read_name = strtok(0, "\t");
        char* tax_field = strtok(0, "\t");
        char* read1_len = strtok(0, "\t");
        char* kmer_str = strtok(0, "\n");

        krp->classified = class_field[0] == 'C' ? true : false;
        krp->read_name = read_name;
        krp->taxid = strtoul(tax_field, (void*)0, 10);
        check_ulong_overflow(krp->taxid);

        krp->read1_len = strtoul(read1_len, (void*)0, 10);
        check_ulong_overflow(krp->read1_len);
        krp->read1_kmers = kmc_fill(krp->read1_kmers, kmer_str);
        return krp;
    }
}

void kraken_print(KrakenRec* krp)
{
    printf("classified: %d\n", krp->classified);
    printf("read_name: %s\n", krp->read_name);
    printf("taxid: %lu\n", krp->taxid);
    printf("read_len1: %lu\n", krp->read1_len);
    printf("read_len2: %lu\n", krp->read2_len);
    printf("read1_kmers size: %zu\n", krp->read1_kmers->size);
    printf("read2_kmers size: %zu\n", krp->read2_kmers->size);
    for (size_t i=0; i < krp->read1_kmers->size; i++){
        printf("%lu:%lu "
                , krp->read1_kmers->taxids[i]
                , krp->read1_kmers->counts[i]
                );
    }
    printf("\n");
    for (size_t i=0; i < krp->read2_kmers->size; i++){
        printf("%lu:%lu "
                , krp->read2_kmers->taxids[i]
                , krp->read2_kmers->counts[i]
                );
    }
    printf("\n");
    return;
}

void kraken_destroy(KrakenRec* krp)
{
    if (krp->read1_kmers)
        kmc_destroy(krp->read1_kmers);
    if (krp->read2_kmers)
        kmc_destroy(krp->read2_kmers);
    free(krp);
    return;
}

FloatHist* fh_create()
{
    FloatHist* fh = malloc(sizeof(*fh));
    fh->size = 1001;
    fh->capacity = 1024;
    fh->fraction_counts = malloc(sizeof(*fh->fraction_counts) * fh->capacity);
    memset(fh->fraction_counts, 0, sizeof(*fh->fraction_counts) * fh->capacity);
    return fh;
}

void fh_destroy(FloatHist* fh)
{
    free(fh->fraction_counts);
    free(fh);
    return;
}

void fh_add(FloatHist* fh, float frac)
{
    int32_t value = (int32_t) round(frac * 1000.0f);
    fh->fraction_counts[value]++;
    /*printf("%f -> %d\n", frac, value);*/
    return;
}

void fh_print(FloatHist* fh)
{
    printf("size: %zu capacity: %zu\n", fh->size, fh->capacity);
    int line = 0;
    printf("%d ", line++);
    for (size_t i=0; i < fh->size; i++){
        printf("%d ", fh->fraction_counts[i]);
        if (!((i+1) % 10)){
            printf("\n%d ",line++);
        }
    }
    return;
}

TaxIdData* txd_create()
{
    TaxIdData* txd = malloc(sizeof(*txd));
    if (txd == 0){
        fprintf(stderr, "Can't allocate TaxidKmerStats\n");
        return txd;
    }
    txd->taxid_size = 0;
    txd->taxid_capacity = TAXID_KMER_DATA_SIZE;
    txd->taxids = malloc(sizeof(*txd->taxids)*TAXID_KMER_DATA_SIZE);

    txd->data = malloc(sizeof(*txd->data)*TAXID_KMER_DATA_SIZE);

    for (int i = 0; i < txd->taxid_capacity; i++){
        txd->data[i] = fh_create();
    }
    if (txd->taxids == 0){
        fprintf(stderr, "Can't allocate TaxidKmerStats\n");
        return txd;
    }
    return txd;
}

void txd_destroy(TaxIdData* txd)
{
    for (int i = 0; i < txd->taxid_capacity; i++){
        fh_destroy(txd->data[i]);
    }
    free(txd->data);
    free(txd->taxids);
    free(txd);
    return;
}

/*//search in O(n)*/
int32_t txd_get_taxa_index(TaxIdData* txd, uint64_t taxid)
{
    int32_t index = -1;
    for (int i=0; i < txd->taxid_size; i++){
        if (txd->taxids[i] == taxid){
            return i;
        }
    }
    return index;
}

TaxIdData* txd_extend(TaxIdData* txd, int32_t new_capacity)
{
    txd->taxid_capacity = new_capacity;
    txd->taxids = realloc(txd->taxids, sizeof(*txd->taxids)*txd->taxid_capacity);
    txd->data = realloc(txd->data, sizeof(*txd->data)*txd->taxid_capacity);
    if (!txd->taxids || !txd->data){
        fprintf(stderr, "Failed to realloc taxids\n");
        return 0;
    } else {
        for (int i = txd->taxid_size; i < txd->taxid_capacity; i++){
            txd->data[i] = fh_create();
        }
        return txd;
    }
}

int32_t txd_add_new_taxa(TaxIdData* txd, uint64_t const taxid)
{
    if (!txd){
        fprintf(stderr, "txd is null\n");
        return -1;
    }
    if (txd->taxid_size+1 == txd->taxid_capacity){
        txd = txd_extend(txd, txd->taxid_capacity + TAXID_KMER_DATA_SIZE);
        if (!txd){
            fprintf(stderr, "Failed to realloc TaxIdData*\n");
            return -1;
        }
    }
    txd->taxids[txd->taxid_size] = taxid;
    return txd->taxid_size++;
}

void txd_add_data(TaxIdData* txd, uint64_t taxid, float frac)
{
    if (frac - 1.0f > 0.00001f){
        fprintf(stderr, "%f is bigger than 1.0f\n", frac);
    }
    int32_t taxid_idx = txd_get_taxa_index(txd, taxid);
    if (taxid_idx == -1){
        taxid_idx = txd_add_new_taxa(txd, taxid);
    }
    fh_add(txd->data[taxid_idx], frac);
    return;
}

void txd_print(TaxIdData* txd)
{

    printf("taxid_size: %d ", txd->taxid_size);
    printf("taxid_capacity: %d\n", txd->taxid_capacity);
    for (int i = 0; i < txd->taxid_size; i++){
        printf("%d taxid: %lu\n",i, txd->taxids[i]);
        fh_print(txd->data[i]);
    }
    printf("\n");
    return;
}

int64_t fh_sum(FloatHist* fh)
{
    int64_t total_kmers = 0;
    for (size_t i = 0; i < fh->size; i++){
        if (total_kmers > INT32_MAX - fh->fraction_counts[i]){
            fprintf(stderr, "Integer overflow at %d\n", __LINE__);
        }
        total_kmers += (int64_t) fh->fraction_counts[i];
    }
    return total_kmers;
}

Quartiles get_quartiles_from_odd(FloatHist* fh, int64_t total_kmers)
{
    Quartiles result = {-1.0f, -1.0f, -1.0f};
    total_kmers += 1;
    float q1_thresh = ((float) (total_kmers / 4) + (total_kmers % 4) / 4.0f);
    float q2_thresh = 2*q1_thresh;
    float q3_thresh = 3*q1_thresh;
    q1_thresh = ceil(q1_thresh);
    q3_thresh = floor(q3_thresh);
    int64_t cumsum = 0;
    for ( size_t i = 0; i < fh->size; i++){
        cumsum += fh->fraction_counts[i];
        if (cumsum >= q1_thresh && result.q1 == -1.0f) {
            result.q1 = i;
        }
        if (cumsum >= q2_thresh && result.q2 == -1.0f) {
            result.q2 = i;
        }
        if (cumsum >= q3_thresh && result.q3 == -1.0f) {
            result.q3 = i;
        }
    }
    return result;
}

Quartiles get_quartiles_from_even(FloatHist* fh, int64_t total_kmers)
{
    Quartiles result = {-1.0f, -1.0f, -1.0f};
    int64_t cumsum = 0;
    total_kmers += 1;

    float q11 = -1.0f;
    float q12 = -1.0f;
    float q21 = -1.0f;
    float q22 = -1.0f;
    float q31 = -1.0f;
    float q32 = -1.0f;
    float thresh_step = (float) (total_kmers / 4) + (total_kmers % 4) / 4.0f;
    float q1_thresh1 = floor(thresh_step);
    float q1_thresh2 = ceil(thresh_step);
    float q2_thresh1 = floor(2.0f*thresh_step);
    float q2_thresh2 = ceil(2.0f*thresh_step);
    float q3_thresh1 = floor(3.0f*thresh_step);
    float q3_thresh2 = ceil(3.0f*thresh_step);
    float q_vals[6] = {q11, q12, q21, q22, q31, q32};
    float q_threshes[6] = {
        q1_thresh1
        , q1_thresh2
        , q2_thresh1
        , q2_thresh2
        , q3_thresh1
        , q3_thresh2
    };
    printf("%f\t%f\t%f\t%f\t%f\t%f\n",
        q1_thresh1
        , q1_thresh2
        , q2_thresh1
        , q2_thresh2
        , q3_thresh1
        , q3_thresh2
        );
    for ( size_t i = 0; i < fh->size; i++){
        cumsum += fh->fraction_counts[i];
        for (size_t j = 0; j < 6; j++){
            if (cumsum >= q_threshes[j] && q_vals[j] == -1.0f) {
                q_vals[j] = i;
                printf("q_vals[%ld] = %ld\n", j, i);
            }
        }
    }
    for (size_t j = 1; j < 6; j++){
        if (q_vals[j] == -1.0f) {
            q_vals[j] = q_vals[j-1];
            printf("q_vals[%ld] = %f\n", j, q_vals[j-1]);
        }

    }


    result.q1 = (q_vals[0] + q_vals[1]) / 2.0f;
    result.q2 = (q_vals[2] + q_vals[3]) / 2.0f;
    result.q3 = (q_vals[4] + q_vals[5]) / 2.0f;
    return result;
}

Quartiles get_nearest_rank_quartiles(FloatHist* fh)
{
    Quartiles result = {-1.0f, -1.0f, -1.0f};
    int64_t total_kmers = fh_sum(fh);

    int64_t q1_thresh = ceil(25*( (float) total_kmers)/100.0f);
    int64_t q2_thresh = ceil(50*( (float) total_kmers)/100.0f);
    int64_t q3_thresh = ceil(75*( (float) total_kmers)/100.0f);

    int64_t cumsum = 0;
    for ( size_t i = 0; i < fh->size; i++){
        cumsum += fh->fraction_counts[i];
        if (cumsum >= q1_thresh && result.q1 == -1.0f) {
            result.q1 = i;
        }
        if (cumsum >= q2_thresh && result.q2 == -1.0f) {
            result.q2 = i;
        }
        if (cumsum >= q3_thresh && result.q3 == -1.0f) {
            result.q3 = i;
        }
    }
    return result;
}

Quartiles get_quartiles(FloatHist* fh)
{
    return get_nearest_rank_quartiles(fh);
    /*Quartiles result = {-1.0f, -1.0f, -1.0f};*/
    /*int64_t total_kmers = fh_sum(fh);*/
    /*result = total_kmers % 2 ? get_quartiles_from_odd(fh, total_kmers) : get_quartiles_from_even(fh, total_kmers);*/
    /*return result;*/
}

float kmer_fraction(KmerCounts const* const kmcs, uint64_t const taxid)
{
    int32_t total_kmers = 0;

    for (size_t i=0; i < kmcs->size; i++){
        total_kmers += kmcs->counts[i];
    }

    float taxid_kmer_fraction = 0.0f;
    for (size_t i=0; i < kmcs->size; i++){
        if (kmcs->taxids[i] == taxid){
            taxid_kmer_fraction += kmcs->counts[i] / (float) total_kmers;
        }
    }

    return taxid_kmer_fraction;
}

float get_avg_kmer_fraction(KrakenRec* krp)
{
    float kmer_frac1 = krp->read1_kmers->size ? kmer_fraction(krp->read1_kmers, krp->taxid) : -1.0f;
    float kmer_frac2 = krp->read2_kmers->size ? kmer_fraction(krp->read2_kmers, krp->taxid) : -1.0f;
    float avg_kmer_frac = kmer_frac1 == -1.0f && kmer_frac2 == -1.0f ? 0.0f :
        (kmer_frac1 != -1.0f && kmer_frac2 == -1.0f) ? kmer_frac1 :
        (kmer_frac2 != -1.0f && kmer_frac1 == -1.0f) ? kmer_frac2 :
        ((kmer_frac1 + kmer_frac2) / 2.0f);
    return avg_kmer_frac;
}

KrakenRec* kraken_adjust_taxonomy(KrakenRec* krp, Taxonomy const* const tx)
{
    if (krp->paired){
        for (size_t i=0; i < krp->read1_kmers->size; i++){
            if (is_a_parent_of_b(krp->taxid, krp->read1_kmers->taxids[i], tx)){
                krp->read1_kmers->taxids[i] = krp->taxid;
            }
        }
        for (size_t i=0; i < krp->read2_kmers->size; i++){
            if (is_a_parent_of_b(krp->taxid, krp->read2_kmers->taxids[i], tx)){
                krp->read2_kmers->taxids[i] = krp->taxid;
            }
        }
    } else {
        for (size_t i=0; i < krp->read1_kmers->size; i++){
            if (is_a_parent_of_b(krp->taxid, krp->read1_kmers->taxids[i], tx)){
                krp->read1_kmers->taxids[i] = krp->taxid;
            }
        }
    }
    return krp;
}
KrakenRec* kraken_adjust_taxonomy_rtl(KrakenRec* krp, Taxonomy const* const tx)
{
    if (krp->paired){
        for (size_t i=0; i < krp->read1_kmers->size; i++){
            if (is_a_parent_of_b(krp->taxid, krp->read1_kmers->taxids[i], tx) || is_a_parent_of_b(krp->read1_kmers->taxids[i], krp->taxid, tx)){
                krp->read1_kmers->taxids[i] = krp->taxid;
            }
        }
        for (size_t i=0; i < krp->read2_kmers->size; i++){
            if (is_a_parent_of_b(krp->taxid, krp->read2_kmers->taxids[i], tx) || is_a_parent_of_b(krp->read2_kmers->taxids[i], krp->taxid, tx)){
                krp->read2_kmers->taxids[i] = krp->taxid;
            }
        }
    } else {
        for (size_t i=0; i < krp->read1_kmers->size; i++){
            if (is_a_parent_of_b(krp->taxid, krp->read1_kmers->taxids[i], tx) || is_a_parent_of_b(krp->read1_kmers->taxids[i], krp->taxid, tx)){
                krp->read1_kmers->taxids[i] = krp->taxid;
            }
        }
    }
    return krp;
}
KmerFractions kmf_calculate(KrakenRec const* const krp)
{
    float kmer_frac1 = -1.0f;
    float kmer_frac2 = -1.0f;
    float avg_kmer_frac = -1.0f;
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
	KmerFractions result = {krp->paired, kmer_frac1, kmer_frac2, avg_kmer_frac};
    return result;
}
