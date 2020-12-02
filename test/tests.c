#include "unity.h"
#include "src/kraken_stats.h"
#include "src/utils.h"

void setUp(void){}
void tearDown(void){}

void test_kmc_create_empty(void)
{
    char str[] = "";
    KmerCounts* kc = kmc_create();
    kc = kmc_fill(kc, str);
    TEST_ASSERT_NOT_NULL(kc);
    TEST_ASSERT_EQUAL(kc->size, 0);
    kmc_destroy(kc);
}

void test_kmc_create_pe_not_empty(void)
{

    char str[] = "1:3 2:3 4:5";
    KmerCounts* kc = kmc_create();
    kc = kmc_fill(kc, str);
    TEST_ASSERT_NOT_NULL(kc);
    kmc_destroy(kc);
}

void test_kmc_create_pe_size_3(void)
{
    char str[] = "1:3 2:3 4:5";
    KmerCounts* kc = kmc_create();
    kc = kmc_fill(kc, str);
    TEST_ASSERT_EQUAL(kc->size, 3);
    kmc_destroy(kc);
}

void test_kmc_create_ambiguous(void)
{
    char str[] = "1:3 2:3 4:5 A:4";
    KmerCounts* kc = kmc_create();
    kc = kmc_fill(kc, str);
    TEST_ASSERT_EQUAL(kc->size, 3);
    kmc_destroy(kc);
}

void test_kmc_create_pe_size_1(void)
{
    char str[] = "4:5";
    KmerCounts* kc = kmc_create();
    kc = kmc_fill(kc, str);
    TEST_ASSERT_EQUAL(kc->size, 1);
    kmc_destroy(kc);
}

void test_kmc_size_3_sum_11(void)
{
    char str[] = "1:3 2:3 4:5";
    KmerCounts* kc = kmc_create();
    kc = kmc_fill(kc, str);
    int sum = 0;
    for (int i = 0; i < kc->size; i++){
        sum += kc->counts[i];
    }
    TEST_ASSERT_EQUAL(sum, 11);
    kmc_destroy(kc);
}

void test_kmc_create_pe_reset_size_0()
{
    char str[] = "1:3 2:3 4:5";
    KmerCounts* kc = kmc_create();
    kc = kmc_fill(kc, str);
    kc = kmc_reset(kc);
    TEST_ASSERT_EQUAL(kc->size, 0);
    kmc_destroy(kc);
}

void test_kraken_empty()
{
    String* str = string_create_from("");
    KrakenRec* krp = kraken_create(true);
    krp = kraken_fill(krp, str);
    /*char empty_str[] = "";*/
    TEST_ASSERT_EQUAL(krp->classified, false);
    TEST_ASSERT_EQUAL_STRING(krp->read_name, "");
    TEST_ASSERT_EQUAL(krp->taxid, 0);
    TEST_ASSERT_EQUAL(krp->read1_len, 0);         
    TEST_ASSERT_EQUAL(krp->read2_len, 0);                                
    TEST_ASSERT_EQUAL(krp->read1_kmers->size, 0);
    TEST_ASSERT_EQUAL(krp->read2_kmers->size, 0);
	string_destroy(str);
    kraken_destroy(krp);
}

void test_kraken_fill_normal()
{
    char s[] = "C\tV100006845L1C001R001118987\t1357\t94|94\t1300:29 1357:26 0:5 |:| 1357:5 0:13 1357:3 0:15 1300:24";
    String* str = string_create_from(s);
    KrakenRec* krp = kraken_create(true);
    krp = kraken_fill(krp, str);
    TEST_ASSERT_EQUAL(krp->classified, true);
    TEST_ASSERT_EQUAL_STRING(krp->read_name, "V100006845L1C001R001118987");
    TEST_ASSERT_EQUAL(krp->taxid, 1357);
    TEST_ASSERT_EQUAL(krp->read1_len, 94);         
    TEST_ASSERT_EQUAL(krp->read2_len, 94);                                
    TEST_ASSERT_EQUAL(krp->read1_kmers->size, 3);
    TEST_ASSERT_EQUAL(krp->read2_kmers->size, 5);
	string_destroy(str);
    kraken_destroy(krp);
}

void test_kraken_fill_read1_empty()
{
    char s[] = "C\tV100006845L1C001R001118987\t1357\t94|94\t |:| 1357:5 0:13 1357:3 0:15 1300:24";
    String* str = string_create_from(s);
    KrakenRec* krp = kraken_create(true);
    krp = kraken_fill(krp, str);
    TEST_ASSERT_EQUAL(krp->classified, true);
    TEST_ASSERT_EQUAL_STRING(krp->read_name, "V100006845L1C001R001118987");
    TEST_ASSERT_EQUAL(krp->taxid, 1357);
    TEST_ASSERT_EQUAL(krp->read1_len, 94);         
    TEST_ASSERT_EQUAL(krp->read2_len, 94);                                
    TEST_ASSERT_EQUAL(krp->read1_kmers->size, 0);
    TEST_ASSERT_EQUAL(krp->read2_kmers->size, 5);
	string_destroy(str);
    kraken_destroy(krp);
}

void test_kraken_fill_read2_empty()
{
    char s[] = "C\tV100006845L1C001R001118987\t1357\t94|94\t1300:29 1357:26 0:5 |:|";
    String* str = string_create_from(s);
    KrakenRec* krp = kraken_create(true);
    krp = kraken_fill(krp, str);
    TEST_ASSERT_EQUAL(krp->classified, true);
    TEST_ASSERT_EQUAL_STRING(krp->read_name, "V100006845L1C001R001118987");
    TEST_ASSERT_EQUAL(krp->taxid, 1357);
    TEST_ASSERT_EQUAL(krp->read1_len, 94);         
    TEST_ASSERT_EQUAL(krp->read2_len, 94);                                
    TEST_ASSERT_EQUAL(krp->read1_kmers->size, 3);
    TEST_ASSERT_EQUAL(krp->read2_kmers->size, 0);
	string_destroy(str);
    kraken_destroy(krp);
}
 
void test_kraken_fill_taxid_overflow()
{
    char s[] = "C\tV100006845L1C001R001118987\t999999999999999999999\t94|94\t1300:29 1357:26 0:5 |:| 1357:5 0:13 1357:3 0:15 1300:24";
    String* str = string_create_from(s);
    KrakenRec* krp = kraken_create(true);
    krp = kraken_fill(krp, str);
    TEST_ASSERT_EQUAL(krp->classified, true);
    TEST_ASSERT_EQUAL_STRING(krp->read_name, "V100006845L1C001R001118987");
    TEST_ASSERT_EQUAL(krp->taxid, 1);
    TEST_ASSERT_EQUAL(krp->read1_len, 94);
    TEST_ASSERT_EQUAL(krp->read2_len, 94);
    TEST_ASSERT_EQUAL(krp->read1_kmers->size, 3);
    TEST_ASSERT_EQUAL(krp->read2_kmers->size, 5);
	string_destroy(str);
    kraken_destroy(krp);
}


void test_kraken_reset()
{
    char s[] = "C\tV100006845L1C001R001118987\t10\t94|94\t1300:29 1357:26 0:5 |:| 1357:5 0:13 1357:3 0:15 1300:24";
    String* str = string_create_from(s);
    KrakenRec* krp = kraken_create(true);
    krp = kraken_fill(krp, str);
    krp = kraken_reset(krp);
    /*char empty_str[] = "";*/
    TEST_ASSERT_EQUAL(krp->classified, false);
    TEST_ASSERT_EQUAL_STRING(krp->read_name, "");
    TEST_ASSERT_EQUAL(krp->taxid, 0);
    TEST_ASSERT_EQUAL(krp->read1_len, 0);
    TEST_ASSERT_EQUAL(krp->read2_len, 0);
    TEST_ASSERT_EQUAL(krp->read1_kmers->size, 0);
    TEST_ASSERT_EQUAL(krp->read2_kmers->size, 0);
	string_destroy(str);
    kraken_destroy(krp);
}

void test_fh_add_1()
{
    FloatHist* fh = fh_create();
    fh_add(fh, 0.001f);
    int32_t sum = 0;
    for (int i=0; i < fh->size; i++){
        sum += fh->fraction_counts[i];
    }
    TEST_ASSERT_EQUAL(sum, 1);
    fh_destroy(fh);
}
void test_fh_add_2()
{
    FloatHist* fh = fh_create();
    fh_add(fh, 0.001f);
    fh_add(fh, 0.0f);
    int32_t sum = 0;
    for (int i=0; i < fh->size; i++){
        sum += fh->fraction_counts[i];
    }
    TEST_ASSERT_EQUAL(sum, 2);
    fh_destroy(fh);
}

void txd_add_taxa_2()
{
    TaxIdData* txd = txd_create();
    int32_t tx_idx = txd_add_new_taxa(txd, 13);
    TEST_ASSERT_EQUAL(tx_idx, 0);
    tx_idx = txd_add_new_taxa(txd, 15);
    TEST_ASSERT_EQUAL(tx_idx, 1);
    txd_destroy(txd);
}

void test_txd_add_data()
{
    TaxIdData* txd = txd_create();
    txd_add_data(txd, 13, 0.3f);
    txd_add_data(txd, 13, 0.9f);
    txd_add_data(txd, 18, 0.1f);
    txd_add_data(txd, 10000, 1.0f);
    TEST_ASSERT_EQUAL(txd->taxid_size, 3);
    txd_destroy(txd);
}

void test_avg_kmer_fraction_pe()
{
    KrakenRec* krp = kraken_create(true);
    char l1[] = "C\tV300023792L1C001R0040001245\t1234\t150|118\t1234:100 0:300 |:| 1234:100 0:200";
    String* str = string_create_from(l1);
    krp = kraken_fill(krp, str);
    float avg_kmer_frac1 = get_avg_kmer_fraction(krp);
    TEST_ASSERT_EQUAL_FLOAT(avg_kmer_frac1, 0.2916667f);
	string_destroy(str);
    kraken_destroy(krp);
}
void test_avg_kmer_fraction_r1_missing_space()
{
    KrakenRec* krp = kraken_create(true);
    char l1[] = "C\tV300023792L1C001R0040001245\t1234\t150|118\t |:| 1234:100 0:200";
    String* str = string_create_from(l1);
    krp = kraken_fill(krp, str);
    float avg_kmer_frac1 = get_avg_kmer_fraction(krp);
    TEST_ASSERT_EQUAL_FLOAT(avg_kmer_frac1, 1.0f/3.0f);
    kraken_destroy(krp);
	string_destroy(str);
}

void test_avg_kmer_fraction_r1_missing_no_space()
{
    KrakenRec* krp = kraken_create(true);
    char l1[] = "C\tV300023792L1C001R0040001245\t1234\t150|118\t|:| 1234:100 0:200";
    String* str = string_create_from(l1);
    krp = kraken_fill(krp, str);
    float avg_kmer_frac1 = get_avg_kmer_fraction(krp);
    TEST_ASSERT_EQUAL_FLOAT(avg_kmer_frac1, 1.0f/3.0f);
	string_destroy(str);
    kraken_destroy(krp);
}

void test_avg_kmer_fraction_r2_missing_space()
{
    KrakenRec* krp = kraken_create(true);
    char l1[] = "C\tV300023792L1C001R0040001245\t1234\t150|118\t1234:100 0:300 |:| ";
    String* str = string_create_from(l1);
    krp = kraken_fill(krp, str);
    float avg_kmer_frac1 = get_avg_kmer_fraction(krp);
    TEST_ASSERT_EQUAL_FLOAT(avg_kmer_frac1, 0.25f);
    kraken_destroy(krp);
	string_destroy(str);
}

void test_avg_kmer_fraction_r2_missing_no_space()
{
    KrakenRec* krp = kraken_create(true);
    char l1[] = "C\tV300023792L1C001R0040001245\t1234\t150|118\t1234:100 0:300 |:|";
    String* str = string_create_from(l1);
    krp = kraken_fill(krp, str);
    float avg_kmer_frac1 = get_avg_kmer_fraction(krp);
    TEST_ASSERT_EQUAL_FLOAT(avg_kmer_frac1, 0.25f);
    kraken_destroy(krp);
	string_destroy(str);
}

void test_quartiles_1_value()
{
    TaxIdData* txd = txd_create();
    txd_add_data(txd, 1, 0.1f);
    Quartiles qs = get_quartiles(txd->data[0]);
    TEST_ASSERT_EQUAL(qs.q1, 100);
    TEST_ASSERT_EQUAL_FLOAT(qs.q2, 100);
    TEST_ASSERT_EQUAL_FLOAT(qs.q3, 100);
    txd_destroy(txd);
}

void test_quartiles_1_value_0()
{
    TaxIdData* txd = txd_create();
    txd_add_data(txd, 1, 0.0f);
    Quartiles qs = get_quartiles(txd->data[0]);
    TEST_ASSERT_EQUAL(qs.q1, 0);
    TEST_ASSERT_EQUAL_FLOAT(qs.q2, 0);
    TEST_ASSERT_EQUAL_FLOAT(qs.q3, 0);
    txd_destroy(txd);
}

void test_quartiles_2_values()
{
    TaxIdData* txd = txd_create();
    txd_add_data(txd, 1, 0.1f);
    txd_add_data(txd, 1, 0.75f);
    Quartiles qs = get_quartiles(txd->data[0]);
    TEST_ASSERT_EQUAL(qs.q1, 100);
    TEST_ASSERT_EQUAL_FLOAT(qs.q2, 100);
    TEST_ASSERT_EQUAL_FLOAT(qs.q3, 750);
    txd_destroy(txd);
}

void test_quartiles_2_values_with_0()
{
    TaxIdData* txd = txd_create();
    txd_add_data(txd, 1, 0.0f);
    txd_add_data(txd, 1, 0.75f);
    Quartiles qs = get_quartiles(txd->data[0]);
    TEST_ASSERT_EQUAL(qs.q1, 0);
    TEST_ASSERT_EQUAL_FLOAT(qs.q2, 0);
    TEST_ASSERT_EQUAL_FLOAT(qs.q3, 750);
    txd_destroy(txd);
}

void test_quartiles_3_values()
{
    TaxIdData* txd = txd_create();
    txd_add_data(txd, 1, 0.1f);
    txd_add_data(txd, 1, 0.2f);
    txd_add_data(txd, 1, 0.3f);
    Quartiles qs = get_quartiles(txd->data[0]);
    TEST_ASSERT_EQUAL(qs.q1, 100);
    TEST_ASSERT_EQUAL_FLOAT(qs.q2, 200);
    TEST_ASSERT_EQUAL_FLOAT(qs.q3, 300);
    txd_destroy(txd);
}

void test_quartiles_3_values_all_0()
{
    TaxIdData* txd = txd_create();
    txd_add_data(txd, 1, 0.0f);
    txd_add_data(txd, 1, 0.0f);
    txd_add_data(txd, 1, 0.0f);
    Quartiles qs = get_quartiles(txd->data[0]);
    TEST_ASSERT_EQUAL(qs.q1, 0);
    TEST_ASSERT_EQUAL_FLOAT(qs.q2, 0);
    TEST_ASSERT_EQUAL_FLOAT(qs.q3, 0);
    txd_destroy(txd);
}

void test_quartiles_1()
{
    TaxIdData* txd = txd_create();
    txd_add_data(txd, 1, 0.1f);
    txd_add_data(txd, 1, 0.2f);
    txd_add_data(txd, 1, 0.3f);
    txd_add_data(txd, 1, 0.4f);
    txd_add_data(txd, 1, 0.5f);
    txd_add_data(txd, 1, 0.6f);
    txd_add_data(txd, 1, 0.7f);
    Quartiles qs = get_quartiles(txd->data[0]);
    TEST_ASSERT_EQUAL(qs.q1, 200);
    TEST_ASSERT_EQUAL_FLOAT(qs.q2, 400);
    TEST_ASSERT_EQUAL_FLOAT(qs.q3, 600);
    txd_destroy(txd);
}
 
void test_quartiles_2()
{
    TaxIdData* txd = txd_create();
    txd_add_data(txd, 1, 0.1f);
    txd_add_data(txd, 1, 0.2f);
    txd_add_data(txd, 1, 0.3f);
    txd_add_data(txd, 1, 0.4f);
    txd_add_data(txd, 1, 0.5f);
    txd_add_data(txd, 1, 0.6f);
    txd_add_data(txd, 1, 0.7f);
    txd_add_data(txd, 1, 0.8f);
    Quartiles qs = get_quartiles(txd->data[0]);
    TEST_ASSERT_EQUAL(qs.q1, 200);
    TEST_ASSERT_EQUAL(qs.q2, 400);
    TEST_ASSERT_EQUAL(qs.q3, 600);
    txd_destroy(txd);
}

void test_quartiles_3()
{
    TaxIdData* txd = txd_create();
    txd_add_data(txd, 1, 0.342);
    txd_add_data(txd, 1, 0.342);
    Quartiles qs = get_quartiles(txd->data[0]);
    TEST_ASSERT_EQUAL(qs.q1, 342);
    TEST_ASSERT_EQUAL_FLOAT(qs.q2, 342);
    TEST_ASSERT_EQUAL_FLOAT(qs.q3, 342);
    txd_destroy(txd);
}


void test_quartiles_4()
{
    TaxIdData* txd = txd_create();
    txd_add_data(txd, 1, 0.03);
    txd_add_data(txd, 1, 0.06);
    txd_add_data(txd, 1, 0.07);
    txd_add_data(txd, 1, 0.08);
    txd_add_data(txd, 1, 0.08);
    txd_add_data(txd, 1, 0.10);
    txd_add_data(txd, 1, 0.13);
    txd_add_data(txd, 1, 0.15);
    txd_add_data(txd, 1, 0.16);
    txd_add_data(txd, 1, 0.20);
    Quartiles qs = get_quartiles(txd->data[0]);
    TEST_ASSERT_EQUAL(qs.q1, 70);
    TEST_ASSERT_EQUAL_FLOAT(qs.q2, 80);
    TEST_ASSERT_EQUAL_FLOAT(qs.q3, 150);
    txd_destroy(txd);
}

void test_quartiles_5()
{
    TaxIdData* txd = txd_create();
    txd_add_data(txd, 1, 0.03);
    txd_add_data(txd, 1, 0.06);
    txd_add_data(txd, 1, 0.07);
    txd_add_data(txd, 1, 0.08);
    txd_add_data(txd, 1, 0.08);
    txd_add_data(txd, 1, 0.09);
    txd_add_data(txd, 1, 0.10);
    txd_add_data(txd, 1, 0.13);
    txd_add_data(txd, 1, 0.15);
    txd_add_data(txd, 1, 0.16);
    txd_add_data(txd, 1, 0.20);
    Quartiles qs = get_quartiles(txd->data[0]);
    TEST_ASSERT_EQUAL(qs.q1, 70);
    TEST_ASSERT_EQUAL_FLOAT(qs.q2, 90);
    TEST_ASSERT_EQUAL_FLOAT(qs.q3, 150);
    txd_destroy(txd);
}


void test_quartiles()
{
    KrakenRec* krp = kraken_create(true);
    TaxIdData* txd = txd_create();
    char l1[] = "C\tV300023792L1C001R0040001245\t1234\t150|118\t1234:100 0:300 |:| 1234:100 0:200";
    char l2[] = "C\tV300023792L1C001R0040001245\t1234\t150|118\t1234:500 A:300 0:500 |:| 1234:200 A:20 0:300";
    char l3[] = "C\tV300023792L1C001R0040001245\t1234\t150|118\t1234:300 0:100 |:| 1234:300 0:100";
    String* str1 = string_create_from(l1);
    String* str2 = string_create_from(l2);
    String* str3 = string_create_from(l3);
    krp = kraken_fill(krp, str1);
    float avg_kmer_frac1 = get_avg_kmer_fraction(krp);
    //0.2916667f
    txd_add_data(txd, krp->taxid, avg_kmer_frac1);
    krp = kraken_reset(krp);

    krp = kraken_fill(krp, str2);
    float avg_kmer_frac2 = get_avg_kmer_fraction(krp);
    //0.45f
    txd_add_data(txd, krp->taxid, avg_kmer_frac2);
    krp = kraken_reset(krp);

    krp = kraken_fill(krp, str3);
    float avg_kmer_frac3 = get_avg_kmer_fraction(krp);
    //0.75f
    txd_add_data(txd, krp->taxid, avg_kmer_frac3);
    krp = kraken_reset(krp);
    
    Quartiles qs = get_quartiles(txd->data[0]);
    TEST_ASSERT_EQUAL(qs.q1, 371);
    TEST_ASSERT_EQUAL_FLOAT(qs.q2, 450);
    TEST_ASSERT_EQUAL_FLOAT(qs.q3, 600);
    kraken_destroy(krp);
	string_destroy(str1);
	string_destroy(str2);
	string_destroy(str3);
}

void test_parents()
{
    Taxonomy* tx = tx_create("test_files/taxo.k2d");
    TEST_ASSERT(is_a_parent_of_b(9605, 9606, tx));
    TEST_ASSERT(is_a_parent_of_b(10088, 10090, tx));
    TEST_ASSERT(is_a_parent_of_b(91347, 562, tx));
    TEST_ASSERT(is_a_parent_of_b(1647988, 239934, tx));
    tx_destroy(tx);
}

void after_string_reset_all_bytes_must_be_0()
{
	String* str = string_create();
	str->str[0] = 'a';
	str->str[1] = 'b';
	str->str[2] = 'c';
	str->str[4094] = 'd';
	str->str[4095] = '\0';
	str->size = 4095;
	string_reset(str);
	bool all_zeros = true;
	for (int i = 0; i < str->capacity; i++){
		all_zeros = all_zeros && (str->str[i] == 0);
	}
	TEST_ASSERT(all_zeros);
}

void read_short_line()
{
	gzFile fh = gzopen("test_files/example.short.1024.out.txt", "rb");
	String* str = string_create();
	parse_line(fh, str);
	TEST_ASSERT(str->size == 1024);
	string_destroy(str);
	gzclose(fh);
}
void read_long_line_bufsize()
{
	gzFile fh = gzopen("test_files/example.long.4096.out.txt", "rb");
	String* str = string_create();
	parse_line(fh, str);
	TEST_ASSERT(str->size == 4096);
	TEST_ASSERT(str->str[0] == 'C');
	TEST_ASSERT(str->str[str->size-1] == '2');
	string_destroy(str);
	gzclose(fh);
}

void read_long_line_bufsize_gz()
{
	gzFile fh = gzopen("test_files/example.long.4096.out.txt.gz", "rb");
	String* str = string_create();
	parse_line(fh, str);
	TEST_ASSERT(str->size == 4096);
	TEST_ASSERT(str->str[0] == 'C');
	TEST_ASSERT(str->str[str->size-1] == '2');
	string_destroy(str);
	gzclose(fh);
}

void read_long_line_bufsize_plus_1()
{
	gzFile fh = gzopen("test_files/example.long.4097.out.txt", "rb");
	String* str = string_create();
	parse_line(fh, str);
	TEST_ASSERT(str->size == 4097);
	TEST_ASSERT(str->str[0] == 'C');
	TEST_ASSERT(str->str[str->size-1] == '4');
	string_destroy(str);
	gzclose(fh);
}

void read_long_line_2x_bufsize()
{
	gzFile fh = gzopen("test_files/example.long.8192.out.txt", "rb");
	String* str = string_create();
	parse_line(fh, str);
	TEST_ASSERT(str->size == 8192);
	TEST_ASSERT(str->str[0] == 'C');
	TEST_ASSERT(str->str[str->size-1] == '2');
	string_destroy(str);
	gzclose(fh);
}

void read_arbitrarily_long_line()
{
	gzFile fh = gzopen("test_files/example.long.15084.out.txt", "rb");
	String* str = string_create();
	parse_line(fh, str);
	TEST_ASSERT(str->size == 15084);
	TEST_ASSERT(str->str[0] == 'C');
	TEST_ASSERT(str->str[str->size-1] == '2');
	string_destroy(str);
	gzclose(fh);
}

int main()
{
    UNITY_BEGIN();
    RUN_TEST(test_kmc_create_empty);
    RUN_TEST(test_kmc_create_pe_not_empty);
    RUN_TEST(test_kmc_create_pe_size_3);
    RUN_TEST(test_kmc_create_pe_size_1);
    RUN_TEST(test_kmc_create_pe_reset_size_0);
    RUN_TEST(test_kmc_size_3_sum_11);
    RUN_TEST(test_kraken_empty);
    RUN_TEST(test_kraken_fill_normal);
    RUN_TEST(test_kraken_fill_read1_empty);
    RUN_TEST(test_kraken_fill_read2_empty);
    RUN_TEST(test_kraken_reset);
	/*RUN_TEST(test_kraken_fill_taxid_overflow);*/
    RUN_TEST(test_fh_add_1);
    RUN_TEST(test_fh_add_2);
    RUN_TEST(txd_add_taxa_2);
    RUN_TEST(test_txd_add_data);
    RUN_TEST(test_avg_kmer_fraction_pe);
    RUN_TEST(test_avg_kmer_fraction_r1_missing_space);
    RUN_TEST(test_avg_kmer_fraction_r1_missing_no_space);
    RUN_TEST(test_avg_kmer_fraction_r2_missing_space);
    RUN_TEST(test_avg_kmer_fraction_r2_missing_no_space);
    RUN_TEST(test_quartiles_1_value);
    RUN_TEST(test_quartiles_1_value_0);
    RUN_TEST(test_quartiles_2_values);
    RUN_TEST(test_quartiles_2_values_with_0);
    RUN_TEST(test_quartiles_3_values);
    RUN_TEST(test_quartiles_3_values_all_0);
    RUN_TEST(test_quartiles_1);
    RUN_TEST(test_quartiles_2);
    RUN_TEST(test_quartiles_3);
    RUN_TEST(test_quartiles_4);
    RUN_TEST(test_quartiles_5);
    RUN_TEST(test_parents);
    RUN_TEST(test_kmc_create_ambiguous);
	RUN_TEST(after_string_reset_all_bytes_must_be_0);
	RUN_TEST(read_short_line);
	RUN_TEST(read_long_line_bufsize);
	RUN_TEST(read_long_line_bufsize_plus_1);
	RUN_TEST(read_long_line_2x_bufsize);
	RUN_TEST(read_arbitrarily_long_line);
	RUN_TEST(read_long_line_bufsize_gz);
    /*RUN_TEST(test_quartiles);*/
    return 0;
}
