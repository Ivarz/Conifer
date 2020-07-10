# Building
```
git clone https://github.com/Ivarz/Conifer && cd Conifer
git submodule update --init --recursive
make
```


# Basic usage
To use this tool you need output file from kraken2 (usually ends with `.out.txt`) and taxonomy database file (`taxo.k2d`).
The following command will calculate confidence score for each classified read. For paired end reads average of the two reads is reported.
```
./conifer -i test_files/example.out.txt -d test_files/taxo.k2d
```

To calculate first three quartiles of the confidence score for each assigned taxonomy use `-s` option:
```
./conifer -s -i test_files/example.out.txt -d test_files/taxo.k2d
```
