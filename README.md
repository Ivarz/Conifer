# Dependencies
[gcc](https://gcc.gnu.org/)
[zlib](https://zlib.net/)
[cmake](https://cmake.org/) if building tests

# Building
```
git clone https://github.com/Ivarz/Conifer && cd Conifer
git submodule update --init --recursive
make
```


# Basic usage
To use this tool you need standard output file from [kraken2](https://github.com/DerrickWood/kraken2) and taxonomy database file (`taxo.k2d`).
The following command will calculate confidence score for each classified read. Note that this kind of output does not include header. For paired end reads confidence score for both reads and the average of the two reads is reported. Only classified reads are reported by default.

```
./conifer -i test_files/example.out.txt -d test_files/taxo.k2d
```

|Kraken standard output| read1 confidence score| read2 confidence score| average confidence score|
|---|---|---|---|
|C V100006960L1C001R001000420 853 100\|100 0:16 853:8 1783272:2 748224:2 1783272:2 168384:5 186801:6 0:2 168384:5 0:18 \|:\| 748224:7 0:2 748224:5 0:21 853:4 748224:7 0:5 748224:3 0:12 |  0.1515 |  0.3939 | 0.2727 |

Use `--rtl` option to obtain RTL scores
```
./conifer --rtl -i test_files/example.out.txt -d test_files/taxo.k2d
```

|Kraken standard output |read1 RTL score| read2 RTL score| average RTL score|
|---|---|---|---|
|C V100006960L1C001R001000420 853 100\|100 0:16 853:8 1783272:2 748224:2 1783272:2 168384:5 186801:6 0:2 168384:5 0:18 \|:\| 748224:7 0:2 748224:5 0:21 853:4 748224:7 0:5 748224:3 0:12 | 0.3636 | 0.6364 | 0.5000

Use `--both_scores` option to obtain confidence and RTL scores simultaneously.
```
./conifer --both_scores -i test_files/example.out.txt -d test_files/taxo.k2d
```
|Kraken standard output| read1 confidence score| read2 confidence score| average confidence score|read1 RTL score| read2 RTL score| average RTL score|
|---|---|---|---|---|---|---|
C V100006960L1C001R001000420 853 100\|100 0:16 853:8 1783272:2 748224:2 1783272:2 168384:5 186801:6 0:2 168384:5 0:18 \|:\| 748224:7 0:2 748224:5 0:21 853:4 748224:7 0:5 748224:3 0:12 |  0.1515 |  0.3939 | 0.2727 | 0.3636 | 0.6364 | 0.5000


```
./conifer -i test_files/example.out.txt -d test_files/taxo.k2d
```

To calculate 25th, 50th and 75th percentiles of the confidence score for each assigned taxonomy use `-s` option.
For paired end reads, average score of each pair is summarized.
For the sake of brevity, only first 5 lines of the summary are shown.

```
./conifer -s -i test_files/example.out.txt -d test_files/taxo.k2d
```
| taxon\_name | taxid | reads | P25 | P50 | P75 |
|---|---|---|---|---|---|
|Faecalibacterium prausnitzii  |  853  |  3  |   0.2200 |  0.2730 | 0.4320|
|Anaerobutyricum hallii | 39488 |  1 |      0.5000 | 0.5000 | 0.5000|
|Lachnospiraceae |186803 | 1 |      0.5000 | 0.5000 | 0.5000|
|Clostridiales |  186802 | 3 |      0.4920 | 0.7200 | 1.0000|

Similar report can be generated for RTL scores:
```
./conifer --rtl -s -i test_files/example.out.txt -d test_files/taxo.k2d
```
| taxon\_name | taxid | reads | P25 | P50 | P75 |
|---|---|---|---|---|---|
|Faecalibacterium prausnitzii  | 853  |  3  |    0.3480 |  0.3480 | 0.4320|
|Anaerobutyricum hallii | 39488 | 1   |   0.5000|  0.5000 | 0.5000|
|Lachnospiraceae | 186803 | 1   |   0.5000| 0.5000|  0.5000|
|Clostridiales   | 186802 | 3   |   0.7200| 1.0000|  1.0000|

and simultaneous reporting of both scores:
```
./conifer --both_scores -s -i test_files/example.out.txt -d test_files/taxo.k2d
```
|taxon\_name  |    taxid |  reads |  P25\_conf |       P50\_conf |       P75\_conf |       P25\_rtl| P50\_rtl| P75\_rtl|
|---|---|---|---|---|---|---|---|---|
|Faecalibacterium prausnitzii |    853   |   3  |      0.2200  | 0.2730  | 0.4320  | 0.3480  | 0.3480  | 0.4320 |
|Anaerobutyricum hallii  | 39488   | 1    |    0.5000  | 0.5000  | 0.5000  | 0.5000  | 0.5000  | 0.5000 |
|Lachnospiraceae  |186803  | 1    |    0.5000  | 0.5000  | 0.5000  | 0.5000  | 0.5000  | 0.5000 |
|Clostridiales   | 186802  | 3     |   0.4920  | 0.7200  | 1.0000  | 0.7200  | 1.0000  | 1.0000 |
